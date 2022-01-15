import sys, os, re, h5py, time
from tqdm import tqdm

import itertools, more_itertools

import numpy as np
import matplotlib.pyplot as plt

import physics_util as pu



### The DT5724B digitizer we will be using has a fixed sampling frequency
### that the user cannot adjust. Thus, this frequency is stored as a 
### global variable within the module, still allowing for changes
### if a different digitizer is used eventually
caen_fsamp = 100.0e6




def _count_generator(reader):
    '''
    Honestly, I'm not even sure how this function works, it's magical. It's
    implemented for efficient counting of the number of lines within a large
    text file, so you don't have to open the file into memory, nor loop over
    the lines as iterables.

    Stolen from: <https://pynative.com/python-count-number-of-lines-in-file/>
    and still exists as of 2021/10/07
    '''
    b = reader(1024 * 1024)
    while b:
        yield b
        b = reader(1024 * 1024)




def count_lines(filename):
    '''
    Using the _count_generator function, search for newline characters in 
    the binary stream that is the file of interest.
    '''
    with open(filename, 'rb') as fp:
        c_generator = _count_generator(fp.raw.read)
        n_lines = sum(buffer.count(b'\n') for buffer in c_generator) + 1

    return n_lines






def fit_waveform_baseline(waveform, pulse_start_ind=0, pulse_end_ind=-1):
    '''
    Function to find the baseline of a pulse from a PMT. Without any of 
    the additional arguments specified, it tries to find the "start" and 
    "stop" points of the pulse, and then exclude the pulse itself from 
    the baseline estimation

    INPUTS

        waveform - a 1D numpy array containing the waveform. Sampling
            information is uncessary for this low-level function

        pulse_start_ind - start of the pulse in units of samples. If 
            the pulse shape/duration is known a priori the user can 
            provide that information. If equal to 0, the function will 
            try to automatically find the start

        pulse_end_ind - end of the pulse in units of samples. If -1, 
            the function will try automatic finding

    OUTPUTS

        baseline - baseline of the pulse in whatever units "waveform"
            input was provided in
    '''

    nsamp = len(waveform)
    if pulse_end_ind == -1:
        pulse_end_ind = nsamp

    inds = np.arange(nsamp)
    bool_inds = (inds < pulse_start_ind) * (inds > pulse_end_ind)

    result = pu.fitting.generate_histogram_and_fit_gaussian(
                waveform[pulse_start_ind:pulse_start_ind], bins=10)

    return result['vals'][1:3]






class WaveformContainer:
    '''
    Class to load .txt files from the Caen digitizer, storing all of the
    waveforms for later integration or peak finding
    '''

    def __init__(self, fname='', n_header_lines=7, no_header=False, \
                 record_length=0, first_n_waveform=0):
        '''
        Initializes a class that can optionally be empty so we don't have to 
        immediately load things if we don't want to. Has the same arguments 
        as the load function, so the docstring looks very similar

        INPUTS

            fname - name of the datafile to load, likely a .txt file

            n_header_lines - lines in the .txt file dedicated to a header and 
                present for EACH event

            no_header - boolean to specify if a header is present or not

            record_length - if no_header==True, then the record length needs
                to be known a priori in order to parse anything at all
        '''

        if not fname:
            self.fname = None

            self.waveform_arr = None
            self.n_waveform = None

        else:
            self.load(fname=fname, n_header_lines=n_header_lines, \
                      no_header=no_header, record_length=record_length, \
                      first_n_waveform=first_n_waveform)



    def save(self):
        '''
        Method to save the class to a pre-defined location for later reloading
        so all of the waveforms don't need to be reprocessed. This is useful
        as the pulse maxima finding algorithm can take a while if the waveforms
        are filtered prior to fitting.
        '''



    def load(self, fname, n_header_lines=7, no_header=False, record_length=0,
             first_n_waveform=0):
        '''
        Function to load a .txt file from the Caen, parsing the distinct events
        intelligently, assuming the text header is part of the file. Without a
        header separating the events, this is probably really difficult without
        knowledge of the digitization window length ahead of time. 

        Given that the Caen digitizer (as it's currently used) loads a configuration
        file for a particular acquisition "session", all events within a single
        datafile will have the same number of samples, and thus we only have to
        determine this number once, then the remainder of the file can be parsed.

        INPUTS

            fname - name of the datafile to load, likely a .txt file

            n_header_lines - lines in the .txt file dedicated to a header and 
                present for EACH event

            no_header - boolean to specify if a header is present or not

            record_length - if no_header==True, then the record length needs
                to be known a priori in order to parse anything at all

            first_n_waveform - option to load only a certain number of waveforms
                from a file, maybe to check some code or something. If 0, the
                function should load everything
        '''

        ### Break if there is not enough information to proceed
        if no_header and not record_length:
            raise ValueError('If no header present, need to provide record length')

        ### Save the filename as a class attribute and quickly check how many lines
        ### are in the whole file
        self.fname = fname

        ### Count the number of lines in the file
        n_lines = count_lines(fname)

        ### Check the header for the first event in the file and try to extract
        ### the waveform length from the header
        if not no_header:

            ### Explicitly extract the first header
            first_header_lines = []
            with open(fname) as file_object:
                for i in range(n_header_lines):
                    first_header_lines.append(file_object.readline())

            ### Look for the term "Record Length", using regular expressions in 
            ### case capitalization is not as expecte
            for line in first_header_lines:
                found_length = re.findall(r'[Rr]ecord\s*[Ll]ength:\s*(\d+)', line)
                if found_length:
                    record_length_found = int(found_length[0])
                    break
                if (i == n_header_lines - 1) and not found_length:
                    raise ValueError("No 'record_length' found in file")
        else:
            ### Define the internal variable for the found length, even if we didn't
            ### find it, since this is the number used in the end
            record_length_found = record_length

        ### If record_length was provided, and a header existed, this checks to make sure
        ### the value found matches the value given
        if record_length:
            if int(record_length) != int(record_length_found):
                raise ValueError("Found a 'record_length' ({:d}) in header and it doesn't" \
                                  .format(record_length_found)
                                  + "match the 'record_length' argument given ({:d}."\
                                  .format(record_length))

        ### Full "line" length of a waveform with its corresponding header
        full_length = record_length_found + n_header_lines

        ### Compute the number of waveforms based on the number of lines and
        ### the known header+record length. Sometimes, the Caen puts an empty
        ### extra line at the end, and the code below allows for this happenstance
        if not first_n_waveform:
            n_waveform = n_lines / full_length
            if (n_lines / full_length).is_integer():
                n_waveform = n_lines / full_length
            elif ((n_lines-1) / full_length).is_integer():
                n_waveform = (n_lines-1) / full_length
            else:
                raise AssertionError("File length is not an integer multiple of the sum " \
                                      + "length: (record_length + n_header_lines).")
        else:
            n_waveform = first_n_waveform

        ### Define the class attribute for the number of waveforms
        self.n_waveform = int(n_waveform)
        self.record_length = record_length_found

        ### Define the output array to which the numbers will be saved, explicitly 
        ### calling a datatype to limit memory usage (14-bit digitizer)
        waveform_arr = np.zeros((self.n_waveform, record_length_found), dtype=np.int16)

        ### Open the file as an object/iterator so we don't have to use memory
        ### to store the list of ASCII strings corresponding to the usual text 
        ### file. A few 100MB file will take a few ~GB stored in RAM as strings
        bad_waveforms = []
        with open(fname, 'r') as file_object:

            ### Chunk the iterator to an iterator of iterators
            waveforms_object = more_itertools.chunked(file_object, full_length)

            for i, waveform in enumerate(waveforms_object):

                ### Little break statement in case we chose to limit the number
                ### of waveforms we wanted to load
                if i >= self.n_waveform:
                    break

                ### Convert the individual waveform iterator into a list, slice
                ### out the header, then save the result to our array
                try:
                    waveform_arr[i] = \
                        np.array(list(waveform)[n_header_lines:], dtype=np.int16)
                except:
                    bad_waveforms.append(i)

        ### Remove the 'bad' waveforms. Sometimes the files have an extra line at the
        ### end or something stupid so often the only 'bad' waveform is a 
        ### non-existent one at the end of the file. This code seems general enough
        ### that maybe it will catch other things (like NaNs or something)
        good_waveforms = np.ones(self.n_waveform) > 0
        for bad_waveform in bad_waveforms:
            print( "Found a 'bad' waveform (index {:d})".format(bad_waveform) )
            print( "   starts at (line {:d})".format(bad_waveform*full_length) )
            print( "   ends at   (line {:d})".format((bad_waveform+1)*full_length) )
            good_waveforms[bad_waveform] = False

        self.waveform_arr = waveform_arr[good_waveforms,:]




    def subtract_baseline(self, fit=False, pulse_start_ind=0, \
                          pulse_end_ind=-1, individual=True):

        inds = np.arange(self.record_length)
        bool_inds = (inds < pulse_start_ind) * (inds > pulse_end_ind)

        if not fit:
            if individual:
                self.baseline_arr = \
                        np.mean(self.waveform_arr[:,bool_inds], axis=1)
            else:
                self.baseline_arr = np.ones(self.n_waveform) * \
                        np.mean(self.waveform_arr[:,bool_inds].flatten())
        
        self.baseline_arr = self.baseline_arr.reshape((self.n_waveform,1))    

        self.waveform_arr_nomean = self.waveform_arr.astype(np.float64) \
                                        - self.baseline_arr



    def find_pulse_maxima_mean(self, sample_window=10, presmooth=False, \
                               presmooth_alpha=0.1):

        if sample_window % 2:
            sample_window += 1

        pulse_maxima = []

        print('Finding Pulse Maxima...')
        for ind in tqdm(range(self.n_waveform)):
            waveform = self.waveform_arr_nomean[ind]
            if presmooth:
                waveform_smooth = pu.filtering.ewma(waveform, presmooth_alpha)
                max_ind = np.argmax(waveform_smooth)
            else:
                max_ind = np.argmax(waveform)

            lower = max_ind - int(0.5*sample_window)
            upper = max_ind + int(0.5*sample_window)

            if lower < 0 or upper > self.record_length:
                print('For pulse ({:d}), max val is near the endpoints'.format(ind))
                pulse_maxima.append(0.0)
            else:
                pulse_maxima.append(\
                        np.mean(waveform[lower:upper]) )

        self.pulse_maxima = np.array(pulse_maxima)




    def find_pulse_maxima_gauss(self, sample_window=30):

        if sample_window % 2:
            sample_window += 1

        pulse_maxima = []

        print('Finding Pulse Maxima...')
        for ind in tqdm(range(self.n_waveform)):
            inds = np.arange(sample_window)
            max_ind = np.argmax(self.waveform_arr_nomean[ind])
            lower = max_ind - int(0.5*sample_window)
            upper = max_ind + int(0.5*sample_window)

            if lower < 0 or upper > self.record_length:
                print('For pulse ({:d}), max val is near the endpoints'.format(ind))
                pulse_maxima.append(0.0)

            else:
                result = pu.fitting.fit_gaussian(inds, \
                            self.waveform_arr_nomean[ind,lower:upper])
                pulse_maxima.append(result['vals'][0])

        self.pulse_maxima = np.array(pulse_maxima)


























