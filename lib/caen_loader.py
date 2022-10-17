import os, sys
from tqdm import tqdm

import dill as pickle

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from matplotlib import cm


### The DT5724B digitizer we will be using has a fixed sampling frequency
### that the user cannot adjust. Thus, this frequency is stored as a 
### global variable within the module, still allowing for changes
### if a different digitizer is used eventually
caen_fsamp = 100.0e6


### It also has a fixed number of bits and input dynamic range
###   CAVEAT: the DT5724 can adjust the dynamic range via attenuation
###           but I'm not entirely sure how to set this up so these
###           are fixed parameters now. Maybe they will be dictionaries
###           that can be adjusted by the user eventually
caen_nbit = 14
caen_range = 2.25
adc_fac = caen_range / (2.0**caen_nbit - 1)


### Some of the functions benefit from a random number generator, so we 
### generate an instance of the default NumPy rng for general use
rng = np.random.default_rng()



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
    the binary stream that is the file of interest. This is useful for the
    .txt files, but we shouldn't be using those. I'll leave the function here
    for now but it might get nix'd.

    INPUTS

        filename - str, name of file with lines to count

    OUTPUTS

        n_lines - int, total number of newline characters + 1, i.e.
            the number of lines in the text file.

    '''
    with open(filename, 'rb') as fp:
        c_generator = _count_generator(fp.raw.read)
        n_lines = sum(buffer.count(b'\n') for buffer in c_generator) + 1

    return n_lines



def get_color_map( n, cmap='plasma', log=False, invert=False):
    '''Gets a map of n colors from cold to hot for use in
       plotting many curves.
       
        INPUTS: 

            n - length of color array to make
            
            cmap - color map for final output

            invert - option to invert

        OUTPUTS: 

            outmap - color map in rgba format
    '''

    n = int(n)
    outmap = []

    if log:
        cNorm = colors.LogNorm(vmin=0, vmax=2*n)
    else:
        cNorm = colors.Normalize(vmin=0, vmax=2*n)

    scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cmap)

    for i in range(n):
        outmap.append( scalarMap.to_rgba(2*i + 1) )

    if invert:
        outmap = outmap[::-1]

    return outmap




def fit_waveform_baseline(waveform, pulse_start_ind=10, pulse_end_ind=-1):
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
    bool_inds = (inds < pulse_start_ind) + (inds > pulse_end_ind)

    print('FITTING BASELINE NOT IMPLEMENTED YET')

    return None






class WaveformContainer:
    '''
    Class to load output files from the CAEN digitizer, storing all of the
    waveforms for later integration or peak finding
    '''

    def __init__(self, fname=None, chunk_size=1000, \
                 filetype='wavedump', **kwargs):
        '''
        Initializes a class that can optionally be empty so we don't have to 
        immediately load things if we don't want to. Has the same arguments 
        as the load function, so the docstring looks very similar

        INPUTS

            fname - str, name of the datafile to load, likely a .dat file

            chunk_size - int, number of waveforms to load into memory 
                simultaneously when processing. This parameter has NOT
                been optimized at all really

            filetype - str, specifies which DAQ framework the file came
                from since the headers and file/data structures are 
                VERY different between the two

            kwargs - various keyword args to pass to the "_preload"
                functions. These have default values where needed and are
                extrememly unlikely to change for things like the
                wavedump files. Thus they're kind of "hidden", but not
                really.
        '''

        if filetype not in ['wavedump', 'compass']:
            assert ValueError("binary file type must be either 'wavedump' "\
                              +"or 'compass', indicating how it was acquired.")

        self.filetype = filetype

        ### Populate the "fname" attribute. Most of the time, the data will
        ### remain on disk, aside from processing/plotting steps
        self.fname = fname
        self.chunk_size = chunk_size

        ### We will count the number of waveforms
        self.n_waveform = None
        
        ### Attributes that get optionally populated depending on 
        ### which methods the user choosed to use
        self.waveform_integrals = None
        self.baseline_arr = None

        if self.filetype == 'wavedump':
            self._preload_wavedump_binary(**kwargs)
            self.waveform_loader = self._load_wavedump_binary

        if self.filetype == 'compass':
            self._preload_compass_binary(**kwargs)
            self.waveform_loader = self._load_compass_binary_waveforms




    def _parse_wavedump_header(self):
        '''
        Parse the first header for the record length, and maybe any other
        parameters if they end up being stored as class attributes
        
        Headers are interpreted based on the UM2019_WaveDump_UserManual_rev18
        released by CAEN. Even though the digitizer manual itself seems to 
        imply the bit patterns in the header are more complicated, when 
        loaded as uint32's, everything just works fine
        '''

        ### Extract the first header from the binary file, making use of the
        ### known header size and datatype. The number of samples in the event
        ### is calculated from the byte_length of the event+header, subtracting
        ### the known header size and noting samples are 14-bit so are packed in
        ### a 16-bit container (2 bytes each). The other header properties
        ### don't seem to be particularly useful
        first_header = np.fromfile(self.fname, dtype='uint32', \
                        offset=0, count=int(0.25*self.header_binary_length))
        byte_length = first_header[0]    ### Header + data
        board_id = first_header[1]       ### Who knows
        pattern = first_header[2]        ### Supposedly a VME board thing
        channel = first_header[3]        ### Digitized channel
        event_counter = first_header[4]  ### Should be '0' for first header

        ### Store the number of samples as a class attribute
        self.record_length = int(0.5*(byte_length - self.header_binary_length))




    def _preload_wavedump_binary(self, header=True, record_length=None, \
                                 n_header_lines=7, header_binary_length=24):
        '''
        Function to "pre-load" a binary file from wavedump. This basically
        parses the header and builds some class attributes that are useful
        to have for some of the later functions

        INPUTS

            record_length - int, if header==False, then the record length 
                needs to be known a priori in order to parse anything at all

            header - bool, specify if a header is present or not

            n_header_lines - int, lines in the .txt file dedicated to a header 
                and present for EACH event

            header_binary_length - int, bytes in the binary .dat file dedicated 
                to the header for EACH event

        OUTPUTS

            NONE - just new class attributes
        '''

        self.header = header

        ### If the file doesn't have a header, a "record_length" argument
        ### has to be passed in order to parse events properly
        if not self.header and self.record_length is None:
            raise ValueError('If no header, need to provide "record_length"')

        ### If there's a header, save the expected size/shape of the header
        ### as class attributes. Again, unlikely that these change or need
        ### to be adjusted from default values.
        if self.header:
            self.n_header_lines = n_header_lines
            self.header_binary_length = header_binary_length
            self._parse_wavedump_header()
        else:
            self.n_header_lines = 0
            self.header_binary_length = 0
            self.record_length = record_length

        ### Compute the total number of waveforms based on the size of the
        ### binary file. This definitely won't work for the ASCII encoded 
        ### files so fuck me on that one. 
        bytes_per_event = self.header_binary_length + 2*self.record_length
        self.n_waveform = int( np.floor(os.path.getsize(self.fname) \
                                          / bytes_per_event) )

        ### Instantiate these arrays as all zeros, so that when chunks
        ### of data are loaded and these arrays are populated, it serves as
        ### a kind of secondary record of what data has been loaded
        self.indices = np.zeros(self.n_waveform)
        self.times = np.zeros(self.n_waveform)

        self.utimes = None




    def _load_wavedump_binary(self, first_index=0, last_index=None, \
                              verbose=True):
        '''
        Function to load a .dat file from the CAEN, assumed to be in binary 
        format. Per the CAEN manual, data are encoded into 2 bytes, even for 
        14-bit digitizers like the one we have. Thus, we can use NumPy's 
        very basic fromfile() function to read the data into an array with 
        dtype=uint16.

        Only loads a small subset of all the waveforms for processing, so
        that they're effectively processed in chunks

        INPUTS

            first_index - int, first event within the file to load, with
                standard python zero-indexing

            last_index - int, last event within the file to load. note 
                the 'last_index' is non-inclusive, such that if
                first_index = 0 and last_index = 1, only the first 
                event will be loaded.

        OUTPUTS

            new_indices - np.ndarray of ints, indices of the returned events

            waveform_arr - np.ndarray of uint16, containing all returned events
                with shape (n_waveform, record_length). Information from the 
                event headers is also stored in arrays, including their unique
                index within the parent file and time from the internal counter
        '''

        ### Figure out the number of waveforms in the current chunk
        if last_index is not None:
            if last_index > self.n_waveform:
                last_index = self.n_waveform
            n_waveform = last_index - first_index
        else:
            n_waveform = self.chunk_size

        ### Figure out the number of bytes to read (maybe store this number
        ### as a class attribute?), but the 'record_length' is necessary by itself
        ### in other functions, as is the header_binary_length
        bytes_per_event = self.header_binary_length + 2*self.record_length
        uint16_to_read = int(0.5 * n_waveform * bytes_per_event )
        assert not bytes_per_event % 2, \
                    "Reading an odd number of bytes. Dun fuq'd up"

        ### Load the 14-bit data from the binary file, noting that CAEN packs
        ### the final result into two bytes, such that we can pretend it's
        ### just a bunch of unsigned 16-bit integers
        raw_data = np.fromfile(self.fname, offset=first_index*bytes_per_event, 
                               count=uint16_to_read, dtype='uint16')

        ### Sometimes, a chunck near the end of the file doesn't have the 
        ### right amount of data, i.e. less than self.chunk_size or less than
        ### last-first, depending on how the function was called. So 
        ### to make sure the array reshaping works on that step, redefine
        ### the quantity n_waveform based on the length of the raw_data 
        ### vector, since np.fromfile() will just read as many bytes are
        ### available if the offset and count parameters are such that 
        ### we're near the end of the file
        if len(raw_data) < uint16_to_read:
            n_waveform = int( len(raw_data) / int(0.5*bytes_per_event) )

        ### Reshape the output so the first axis indexes the waveforms, and
        ### the second axis indexes time
        waveform_arr = raw_data.reshape(n_waveform, int(0.5*bytes_per_event))

        if self.header:
            ### Parse the header and extract the event indices and event times,
            ### to be incorporated into the indexing and time-keeping class
            ### attributes defined to be nominally empty in the init() statement
            new_indices = np.frombuffer(waveform_arr[:,8:10].flatten().tobytes(), \
                                        dtype='uint32')
            new_times = np.frombuffer(waveform_arr[:,10:12].flatten().tobytes(), \
                                      dtype='uint32')

            self.indices[new_indices] = new_indices

            ### 'ADJUSTING' THE TIME STAMP BASED ON MANUAL
            ### Shift left by one bit, then rignt by one bit since the digitizer
            ### has a 31-bit counter for time stamps, with an overflow bit that
            ### we will ignore (it's 0 for the first count, then continuously 1
            ### until the next acquisition session, i.e. useless). This left 
            ### then right shift effectively clears the MSB. Also, convert from
            ### 10-ns clock cycles to actual seconds
            self.times[new_indices] = \
                10.0e-9 * np.right_shift(np.left_shift(new_times, 1), 1)

        return new_indices, waveform_arr[:,int(0.5*self.header_binary_length):]





    def _parse_compass_header(self):
        '''
        Parse the initial file header from the CoMPASS binary to see what 
        kind of things are in the file. 
        '''

        with open(self.fname, 'rb') as file:
            first_byte = int.from_bytes(file.read(1), byteorder='little')

        out = {}
        out['energy'] = (first_byte & 0b00000001) > 0
        out['keV/MeV'] = (first_byte & 0b00000010) > 0
        out['energy_short'] = (first_byte & 0b00000100) > 0
        out['waveforms'] = (first_byte & 0b00001000) > 0

        return out




    def _preload_compass_binary(self):
        '''
        Function to "pre-load" a binary file from wavedump. This basically
        parses the header and builds some class attributes that are useful
        to have for some of the later functions
        '''

        self.compass_binary_contents = self._parse_compass_header()
        self.header_binary_length = 16

        if self.compass_binary_contents['energy']:
            self.header_binary_length += 2
        
        if self.compass_binary_contents['keV/MeV']:
            self.header_binary_length += 8
        
        if self.compass_binary_contents['energy_short']:
            self.header_binary_length += 2
        
        if self.compass_binary_contents['waveforms']:
            self.header_binary_length += 5

        with open(self.fname, 'rb') as file:
            ### Skip the overall file header (2 bytes) and the first event
            ### header of the given length. Go back 4 bytes to get the 
            ### number of waveforms (hence the pedantic maths below).
            ### The number of samples is stored as a 32-bit integer.
            file.seek(2 + self.header_binary_length - 4)
            self.record_length = int.from_bytes(file.read(4), byteorder='little', \
                                                signed=False)

        ### Compute the total number of waveforms based on the size of the
        ### binary file. This definitely won't work for the ASCII encoded 
        ### files so fuck me on that one. 
        bytes_per_event = self.header_binary_length + 2*self.record_length
        self.n_waveform = int( np.floor( (os.path.getsize(self.fname) - 2) \
                                          / bytes_per_event ) )

        ### Instantiate these arrays as all zeros, so that when chunks
        ### of data are loaded and these arrays are populated, it serves as
        ### a kind of secondary record of what data has been loaded
        self.indices = np.zeros(self.n_waveform)
        self.times = np.zeros(self.n_waveform)
        self.waveform_types = np.zeros(self.n_waveform)

        self.utimes = None





    def _load_compass_binary_waveforms(\
            self, first_index=0, last_index=None, verbose=True):
        '''
        Function to load a .dat file from the CAEN, assumed to be in binary 
        format. Per the CAEN manual, data are encoded into 2 bytes, even for 
        14-bit digitizers like the one we have. Thus, we can use NumPy's 
        very basic fromfile() function to read the data into an array with 
        dtype=uint16.

        Only loads a small subset of all the waveforms for processing, so
        that they're effectively processed in chunks

        INPUTS

            first_index - int, first event within the file to load, with
                standard python zero-indexing

            last_index - int, last event within the file to load. note 
                the 'last_index' is non-inclusive, such that if
                first_index = 0 and last_index = 1, only the first 
                event will be loaded.

        OUTPUTS

            new_indices - np.ndarray of ints, indices of the returned events

            waveform_arr - np.ndarray of uint16, containing all returned events
                with shape (n_waveform, record_length). Information from the 
                event headers is also stored in arrays, including their unique
                index within the parent file and time from the internal counter
        '''

        ### Figure out the number of waveforms in the current chunk
        if last_index is not None:
            if last_index > self.n_waveform:
                last_index = self.n_waveform
            n_waveform = last_index - first_index
        else:
            n_waveform = self.chunk_size

        ### Figure out the number of bytes to read (maybe store this number
        ### as a class attribute?), but the 'record_length' is necessary by itself
        ### in other functions, as is the header_binary_length
        bytes_per_event = self.header_binary_length + 2*self.record_length
        bytes_to_read = int(n_waveform * bytes_per_event)

        ### Load the raw data with headers from the binary file as uint8s
        ### so each element is a single byte
        raw_data = np.fromfile(self.fname, offset=first_index*bytes_per_event + 2, 
                               count=bytes_to_read, dtype='uint8')

        ### Sometimes, a chunck near the end of the file doesn't have the 
        ### right amount of data, i.e. less than self.chunk_size or less than
        ### last-first, depending on how the function was called. So 
        ### to make sure the array reshaping works on that step, redefine
        ### the quantity n_waveform based on the length of the raw_data 
        ### vector, since np.fromfile() will just read as many bytes are
        ### available if the offset and count parameters are such that 
        ### we're near the end of the file
        if len(raw_data) < bytes_to_read:
            n_waveform = int( len(raw_data) / int(bytes_per_event) )

        ### For the wavedump binaries, the file actually tells you the 
        ### indices of the recorded events, but it's always just counting
        ### up the natural numbers. For the CoMPASS binaries, the
        ### index is implicit from the location within the file
        new_indices = np.arange(n_waveform, dtype=int) + first_index
        self.indices[new_indices] = new_indices

        ### Reshape the output so the first axis indexes the waveforms
        data_arr = raw_data.reshape(n_waveform, bytes_per_event)

        ### Extract the waveform time-series from the raw bytearray, convert
        ### it all to 16-bit unsigned integers and the reshape the flat
        ### array into a 2D array with the first axis again indexing 
        ### waveforms and the second indexing time
        waveform_data = \
            np.frombuffer(data_arr[:,self.header_binary_length:]\
                                .flatten().tobytes(), dtype='int16')
        waveform_arr = waveform_data.reshape(n_waveform, self.record_length)

        ### Parse the header and extract the event times, to be incorporated 
        ### into the time-keeping class attributes defined to be nominally empty 
        ### in the "_preload" functions
        new_times = np.frombuffer(data_arr[:,4:12].flatten().tobytes(), \
                                  dtype='uint64')
        waveform_types = np.frombuffer(data_arr[:,20].tobytes(), dtype='uint8')
        self.waveform_types[new_indices] = waveform_types      

        ### 'ADJUSTING' THE TIME STAMP BASED ON MANUAL
        ### Shift left by one bit, the rignt by one bit since the digitizer
        ### has a 31-bit counter for time stamps, with an overflow bit that
        ### we will ignore (it's 0 for the first count, then continuously 1
        ### until the next acquisition session, i.e. useless). This left 
        ### then right shift effectively clears the MSB. Also, convert from
        ### 10-ns clock cycles to actual seconds
        self.times[new_indices] = 1.0e-12 * new_times

        return new_indices, waveform_arr





    def _load_compass_binary_energies(self, plot=False):
        '''
        Function to load a .dat file from the CAEN, assumed to be in binary 
        format. Per the CAEN manual, data are encoded into 2 bytes, even for 
        14-bit digitizers like the one we have. Thus, we can use NumPy's 
        very basic fromfile() function to read the data into an array with 
        dtype=uint16.

        Only loads a small subset of all the waveforms for processing, so
        that they're effectively processed in chunks

        INPUTS

            plot - boolean, flag for quick debug plotting

        OUTPUTS

        '''

        ### Figure out the number of bytes to read (maybe store this number
        ### as a class attribute?), but the 'record_length' is necessary by itself
        ### in other functions, as is the header_binary_length
        bytes_per_event = self.header_binary_length + 2*self.record_length

        if self.waveform_integrals is None:
            self.waveform_integrals = np.zeros(self.n_waveform)

        with open(self.fname, 'rb') as file:

            ### For each event, skip the overall file header (2 bytes),
            ### advance to the event's position within the file, read the
            ### 16-bit value encoding the energy (bytes n+12 and n+13), 
            ### and keep looping
            for i in range(self.n_waveform):
                index = 2 + i*(bytes_per_event) + 12

                file.seek(index, 0)
                self.waveform_integrals[i] = \
                    int.from_bytes(file.read(2), \
                                   byteorder='little', \
                                   signed=False)

        if plot:
            bins = np.arange(np.max(self.waveform_integrals)+1)
            plt.hist(self.waveform_integrals, bins=bins)
            plt.tight_layout()
            plt.show()





    def _unwrap_times(self, counter_nbit=31, counter_freq=100.0e6):
        '''
        Takes the self.times() array containing the "trigger time tags"
        from each event, and unwraps them, as the internal counter on the
        DT5724 digitizer only counts up to ~21 seconds before the 
        31-bit counter runs out and rolls over.

        There is supposedly an Extended Trigger Time Tag (ETTT) option in the
        software, where some header bits get dedicated to an extra 14 bits of 
        counter, effectively extending the full time to many hours. Haven't
        figured out how to use this feature though.

        Adjusts the self.times class attribute in-place

        INPUTS

            counter_nbit - int, number of bits in the counter, to help
                determine the max time the counter can get up to

            counter_freq - float, frequency of the counter, in Hz

        '''

        ### Determine the total length of one full counter sequence
        counter_time = 1.0 / counter_freq
        total_time = counter_time * 2**counter_nbit

        ### With a simple finite difference method, determine where the timer
        ### wraps over, since it's a really sharp feature. The 'wraps' array
        ### that's built at the end of this are just the indices of the 
        ### self.times array where where the timer loops 
        rate = np.diff(self.times)
        wraps = np.arange(len(rate))[np.where(np.abs(rate) > 0.5*total_time)]

        ### Define a new array to make sure we don't get into any weird
        ### python referencing business. I honestly don't know if I need to 
        ### do this, but when I tried to do in-place modifcations of the 
        ### self.time class attribute, it wasn't working properly
        new_times = self.times.copy()

        ### For each time the timer wraps, add the full duration of the
        ### counter to every point following
        for wrap in wraps:
            new_times[wrap+1:] += total_time

        self.utimes = new_times




    def plot_waveform(self, index, baseline=True, tscale='us', \
                      amp_scale='bits', fig=None, ax=None, \
                      show=True, color=None, **kwargs):
        '''
        Plots a specific waveform, based on the static index assigned to it
        in the raw data file, and preserved in the self.indices class attribute

        INPUTS

            index - int, waveform index to plot

            baseline - boolean, whether to include the baseline when plotting.
                If False, the pre-pulse data should be mean-zero

            tscale - str, scale for the time axis. If equal to 'us', 'u', or 
                'U', the waveform in plotted in microseconds. If equal to 'ns',
                'n' or 'N', the waveform is shown in nanoseconds

            amp_scale - str, scale for the vertical axis. If equal to 'bits',
                'b', or 'B', plot the y-axis in bits. If equal to 'volts', 'v'
                or 'V', plot in volts with the hard-coded conversion factor

            show - boolean, whether to show the figure or not

        OUTPUTS

            fig - matplotlib.pyplot.fig, handle for the figure object.
                Useful to overlay stuff on default plots

            ax - matplotlib.pyplot.axes, handle for the axes object 
                containing all the plotted things
        '''

        if fig is None or ax is None:
            fig, ax = plt.subplots(1,1)

        if color is None:
            color = 'C0'

        ### Make some labels and factors for nice plotting
        if tscale == 'us' or tscale == 'u' or tscale =='U':
            tscale_fac = 1e6
            xlabel = 'Time [$\\mu$s]'
        elif tscale == 'ns' or tscale == 'n' or tscale == 'N':
            tscale_fac = 1e9
            xlabel = 'Time [ns]'
     
        if amp_scale == 'bits' or amp_scale == 'b' or amp_scale == 'B':
            amp_scale_fac = 1.0
            ylabel = 'PMT signal [ADC Counts]'  
        elif amp_scale == 'volts'or amp_scale == 'v' or amp_scale == 'V':
            amp_scale_fac = adc_fac
            ylabel = 'PMT signal [V]'

        ### Load the desired waveform with the binary loading function
        binary_index_arr, raw_waveform_arr = \
            self.waveform_loader(first_index=index, last_index=index+1, 
                                 verbose=False)

        ### Extract the desired event from the naive array that comes from
        ### the binary loading function
        binary_index = binary_index_arr[0]
        raw_waveform = raw_waveform_arr[0]

        ### Build a time vector for plotting
        tvec = np.arange(len(raw_waveform)) * (1.0 / caen_fsamp)

        ### Sass the user if they try to do something bad
        if (not baseline) and (self.baseline_arr is None):
            print('Compute the baseline before trying to plot without it...')
            baseline = True
        waveform = raw_waveform.astype(np.float64).copy()

        ### Optionally subtract the baseline. Also, plot a horizontal line
        ### either at y=0 if the baseline is subtracted, or at the level of
        ### the baseline, for ease of viewing
        if self.baseline_arr is not None:
            baseline_val = self.baseline_arr[index]
            if baseline:
                ax.axhline(baseline_val*amp_scale_fac, color='k', ls='--', \
                           alpha=0.5, zorder=3)
            else:
                waveform -= baseline_val
                ax.axhline(0.0, color='k', ls='--', alpha=0.5, zorder=3)

        ### Plot the waveform
        ax.plot(tvec*tscale_fac, waveform*amp_scale_fac, color=color)

        ### Label things
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        fig.tight_layout()

        if show:
            plt.show()

        ### Returns the figure and axes objects if desired so multiple 
        ### spectra can be plotted over one another
        return fig, ax






    def compute_baseline(self, fit=False, plot=False, pulse_start_ind=10, \
                          pulse_end_ind=-1, individual=True, \
                          bad_threshold_sigma=3.5, **kwargs):
        '''
        Method to compute the pulse baseline, by taking a simple average
        of the waveform near the beginning of the digitized event, and 
        optionally points near the end as well.

        INPUTS

            fit - boolean, to fit a Gaussian instead of taking a simple average
                although this has yet to be implemented. default: False

            plot - boolean, to plot the computed baseline, highlighting events
                where the baseline is far from the mean value of the baseline

            pulse_start_ind - int, lower index (in units of samples) of the 
                actual pulse within the waveform. Baseline is computed from 
                points up to this index. default: 10

            pulse_end_ind - int, upper index of the pulse within the waveform. 
                Can be used to optionally include points from the end in the
                baseline computation. default: -1

            individual - boolean, selects whether to compute the baseline for
                each pulse individually, or compute a common baseline

            bad_threshold_sigma - float, number of standard deviations beyond
                which an individual baseline is considered 'too far' from the
                'normal' behavior of the baseline

            **kwargs - passed to the 'self.plot_baseline()' function

        '''

        print('Computing waveform baselines...')

        if fit:
            raise NotImplementedError(\
                     'Fitting the baseline is not yet an option')

        ### Make an array of indices [0,...,record_length] for easy subselection
        ### of different parts of teh waveform
        inds = np.arange(self.record_length)
        if pulse_end_ind == -1:
            pulse_end_ind = self.record_length
            
        ### Determine which points to use in the baseline estimation
        bool_inds = (inds < pulse_start_ind) + (inds > pulse_end_ind)

        ### Allocate the class attribute that will hold the values of the 
        ### baseline from each event
        self.baseline_arr = np.zeros(self.n_waveform)

        ### Loop over all the events, loading them in chunks to reduce
        ### RAM pressure
        for chunk in range( int(self.n_waveform/self.chunk_size) + 1 ):

            ### Load the chunk
            chunk_index = chunk * self.chunk_size
            indices, waveforms = \
                self.waveform_loader(first_index=chunk_index, \
                                     last_index=chunk_index+self.chunk_size)

            ### Compute the baselines with a simple vector operation
            if not fit:
                if individual:
                    new_baselines = np.mean(waveforms[:,bool_inds], axis=-1)
                else:
                    new_baselines = np.ones(len(indices)) * \
                                        np.mean(waforms[:,bool_inds].flatten())
                self.baseline_arr[indices] = new_baselines

        ### Pre-select 'bad' baselines, if they're too far away from the 
        ### mean value. Usually these are weirdly triggered events, or events
        ### with double pulses or something like that. The event data is still
        ### preserved, but can be excluded by using this boolean indexing array
        mean_baseline = np.mean(self.baseline_arr)
        std_baseline = np.std(self.baseline_arr)
        self.bad_baselines = \
            np.abs(self.baseline_arr - mean_baseline) \
                > bad_threshold_sigma*std_baseline

        ### The event time counter is a 31-bit counter operating at 100 MHz,
        ### so it's max value is ~21 seconds. This internal class method
        ### unwraps the overflowed counter values in some event headers, and 
        ### adjusts the self.times class attribute in place
        self._unwrap_times()

        ### Plot if desired
        if plot:
            self.plot_baseline(**kwargs)




    def plot_baseline(self, tscale='s', amp_scale='bits', \
                      mean_subtract=False, show=True):
        '''
        Plots the computed baseline as a function of time. Mostly a debugging
        tool, or to validate whether data is behaving strangely.

        INPUTS

            tscale - str, scale for the time axis. If equal to 'us', 'u', or 
                'U', the waveform in plotted in microseconds. If equal to 'ns',
                'n' or 'N', the waveform is shown in nanoseconds

            amp_scale - str, scale for the vertical axis. If equal to 'bits',
                'b', or 'B', plot the y-axis in bits. If equal to 'volts', 'v'
                or 'V', plot in volts with the hard-coded conversion factor

            mean_subtract - boolean, whether to show the baseline centered 
                around 0 (maybe a useful thing?), or leave it at the
                actual value computed

            show - boolean, whether to show the figure or not

        OUTPUTS

            fig - matplotlib.pyplot.fig, handle for the figure object.
                Useful to overlay stuff on default plots

            ax - matplotlib.pyplot.axes, handle for the axes object 
                containing all the plotted things
        '''

        if self.baseline_arr is None:
            raise ValueError("Need to compute the baseline before plotting " \
                             + "which may require special arguments that I " \
                             + "can't predict for you.")

        ### Make some labels and factors for nice plotting
        if tscale == 's':
            tscale_fac = 1.0
            xlabel = 'Time [s]'
        elif tscale == 'hr' or tscale == 'h':
            tscale_fac = 1.0 / 3600
            xlabel = 'Time [hr]'

        if amp_scale == 'bits' or amp_scale == 'b' or amp_scale == 'B':
            amp_scale_fac = 1.0
            ylabel = 'Signal Baseline [ADC Counts]'
        elif amp_scale == 'volts'or amp_scale == 'v' or amp_scale == 'V':
            amp_scale_fac = adc_fac
            ylabel = 'Signal Baseline [V]'

        ### Define the figure and axes handles
        fig, ax = plt.subplots(1,1)

        ### Compute the mean and standard deviation of the baseline values for
        ### all events in this file
        good_inds = np.logical_not(self.bad_baselines)
        mean = np.mean(self.baseline_arr[good_inds])
        std = np.std(self.baseline_arr[good_inds])

        ### Determine any plotting offset that might be required, as well as
        ### plot a horizontal line, to assist in visualization, at the value 
        ### of the baseline, or at 0 if the mean has been subtracted
        offset = 0.0
        if mean_subtract:
            offset = -1.0 * mean
            ax.axhline(0.0, color='k', ls='--', alpha=0.5, zorder=3)
        else:
            ax.axhline(mean*amp_scale_fac, color='k', ls='--', \
                       alpha=0.5, zorder=3)
            ax.set_title(f'Mean = {mean*amp_scale_fac:0.6g},  '
                         + f'Std. dev. = {std*amp_scale_fac:0.2g}', fontsize=14)

        ### Plot the stuff
        ax.plot(self.utimes*tscale_fac, \
                (self.baseline_arr+offset)*amp_scale_fac, \
                color='r', alpha=0.5, ls=':', zorder=5)
        ax.plot(self.utimes[good_inds]*tscale_fac, \
                (self.baseline_arr[good_inds]+offset)*amp_scale_fac, \
                color='C0', zorder=5)

        ### Make some labels
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        fig.tight_layout()

        if show:
            plt.show()

        return fig, ax






    def integrate_waveforms(self, adaptive_window=True, asymmetry=0.5, \
                            integration_start=0.0e-6, \
                            integration_window=500.0e-9, \
                            polarity=-1, plot=False, \
                            plot_index=None, **kwargs):
        '''
        Method to integrate the pulse area, via some simple summing operations
        on arrays. Has some options to setup integration windows that are 
        static (assuming good triggering) or dynamic based on the pulse 
        maxima.

        Defines a class-attribute and doesn't immediately return anything
        to the user

        INPUTS

            adaptive_window - boolean, whether to adjust the left and right 
                edges of the window around the most extreme value. This
                argument is default True since a static timing window 
                seems generically dumb

            asymmetry - float [0,1], asymmetry of the adaptive window. A lower
                number means there are fewer time-series points included in the
                integral estimation from BEFORE the extremum value, e.g. a 
                value of 0.3 would mean: 
                    t_start = t_extremum - 0.3*integration_window
                      t_end = t_extremum + (1 - 0.3)*integration_window

            integration_start - float, start time of the integration window

            integration_window - float, length of time to integrate the 
                pulse area

            polarity - +1 or -1, polarity of the pulse. Used to decide on
                either using argmax or argmin to find the extremum

            plot - boolean, whether to plot an example of the integration
                window, useful for debugging

            plot_index - int, index of waveform to plot if plotting is 
                desired. Will generate a random number if none given

            **kwargs

        '''

        if plot:
            if plot_index is None:
                plot_index = rng.integers(low=0, high=self.n_waveform)

        print('Integrating waveforms...')

        ### Time vector for establishing an integration window
        tvec = np.arange(self.record_length) * (1.0 / caen_fsamp)

        ### Load the data in chunks to keep down RAM usage
        self.waveform_integrals = np.zeros(self.n_waveform)
        for chunk in range( int(self.n_waveform/self.chunk_size) + 1 ):

            ### Load only the waveforms in this particular chunk
            chunk_index = chunk * self.chunk_size
            indices, waveforms = \
                self.waveform_loader(first_index=chunk_index, \
                                     last_index=chunk_index+self.chunk_size)

            ### Define this quantity for easy array building. Nominally it should
            ### be the same as self.chunk_size, except for the last chunk.
            ### This just preserves generality regardless
            n_waveform = waveforms.shape[0]

            if adaptive_window:

                ### Build the indexing array and a tiled time array
                full_inds = np.zeros((n_waveform, self.record_length), \
                                     dtype=bool)
                full_times = np.tile(tvec,(n_waveform,1))

                ### Find the pulse extrema, dependent on the polarity
                if polarity > 0:
                    ext_time_inds = np.argmax(waveforms, axis=-1)
                elif polarity < 0:
                    ext_time_inds = np.argmin(waveforms, axis=-1)

                ### Build the integration windows based on the extrema times,
                ### in a vectorized fashion
                ext_times = tvec[ext_time_inds]
                lower_times = (ext_times - asymmetry*integration_window)\
                                    .reshape((n_waveform,1))
                upper_times = (ext_times + (1.0 - asymmetry)*integration_window)\
                                    .reshape((n_waveform,1))

                ### Get an indexing array for the full array of waveforms.
                ### Allows for easy summation at the end of the day
                full_inds = (full_times > lower_times) * (full_times < upper_times)

            else:

                ### If we can reliably integrate the same static time window for
                ### every pulse, we can use the derpy arguments to do just that
                inds = (tvec > integration_start) \
                            * (tvec < integration_start + integration_window)
                full_inds = np.tile(inds,(n_waveform,1))

            ### If we want to plot the integration, it helps to remember
            ### the specific integration window we used for the pulse we want to
            ### plot, since it can be adaptive
            if plot:
                if np.abs(chunk_index - plot_index) < self.chunk_size \
                        and plot_index > chunk_index:
                    plot_window = full_inds[plot_index-chunk_index]
                    plot_ext_time = ext_times[plot_index-chunk_index]

            ### Raises an AttributeError if the baseline hasn't been computed.
            ### One would probably never want to include the baseline in the 
            ### pulse integral estimate
            if self.baseline_arr is None:
                raise AttributeError("You probably want to subtract the baseline "\
                                     + "before you integrate the pulses...")

            ### Up-convert the waveform datatype for naive math with the 
            ### baseline array (since that's necessarily a float)
            waveforms64 = waveforms.astype(np.float64)
            baselines = self.baseline_arr[indices].reshape((n_waveform,1))

            ### Compute the pulse integrals as a simple sum over all points
            ### within the integration window, using the full indexing array built
            ### earlier as a mask. Given the non-uniform way the adaptive window
            ### indexing happens, if you try to use the indices directly, NumPy
            ### is sad and flattens the array.
            self.waveform_integrals[indices] = polarity * \
                np.sum( (waveforms64 - baselines)*full_inds, axis=-1)

        ### Plot that baby if desired!
        if plot:

            ### Use the built in plotting method to generate a figure and axes
            ### object with the waveform to be plotted, but don't show it yet
            ### so we can add some things
            fig, ax = self.plot_waveform(plot_index, tscale='u', \
                                         show=False, **kwargs)

            ### Get the limits of the nice waveform plot so after we muck things
            ### up, we can set it right again
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()

            ### Shade the integration window
            ax.fill_between(tvec*1e6, 0, 1, where=plot_window, \
                            color='k', alpha=0.35, \
                            transform=ax.get_xaxis_transform())

            ### Add a red vertical line at the location of the extremum value
            ax.axvline(plot_ext_time*1e6, color='r', alpha=0.5, ls=':')

            ### Put a little derpy label with the pulse integral
            ax.text(0.7, 0.5, f'Integral = {self.waveform_integrals[plot_index]:0.3g}', \
                    ha='center', va='center', transform=ax.transAxes)

            ### Fix the axes limits and run a tight_layout()
            ax.set_xlim(*xlim)
            ax.set_ylim(*ylim)
            fig.tight_layout()
            plt.show()

        return None




    def find_pulse_extrema(self, npoints=3, polarity=-1):
        '''
        INPUTS

            npoints - int, number of points to average when computing the
                value of the extremum

            polarity - +1 or -1, polarity of the pulse. Used to decide on
                either using argmax or argmin to find the extremum

            plot - boolean, whether to plot an example of the integration
                window, useful for debugging

            plot_index - int, index of waveform to plot if plotting is 
                desired

            **kwargs

        '''

        print('Finding Pulse Extrema...')

        ### Load the data in chunks to keep down RAM usage
        self.waveform_extrema = np.zeros(self.n_waveform)
        for chunk in range( int(self.n_waveform/self.chunk_size) + 1 ):

            ### Load only the waveforms in this particular chunk
            chunk_index = chunk * self.chunk_size
            indices, waveforms = \
                self.waveform_loader(first_index=chunk_index, \
                                     last_index=chunk_index+self.chunk_size)

            ### Define this quantity for easy array building. Nominally it should
            ### be the same as self.chunk_size, except for the last chunk.
            ### This just preserves generality regardless
            n_waveform = waveforms.shape[0]

            inds = np.arange(n_waveform)

            ### Find the pulse extrema, dependent on the polarity
            if polarity > 0:
                ext_time_inds = np.argmax(waveforms, axis=-1)
            elif polarity < 0:
                ext_time_inds = np.argmin(waveforms, axis=-1)

            ### Get the extremum values, and the baseline
            ext_vals = waveforms[inds,ext_time_inds].astype(np.float64)

            baselines = self.baseline_arr[indices]

            ### Add the extrema values, minus the baseline, to the array
            ### taking an explicit absolute value (for sensible spectra
            ### plotting), assuming the user is keeping track of the
            ### pulse polarity
            self.waveform_extrema[indices] = \
                        polarity * (ext_vals - baselines)


        return None






    def plot_pulse_spectra(self, nbin=100, spectra_type='int', hist_range=None, \
                           log_scale=True, filled=True, exclude_bad_baseline=True, \
                           color=None, fig=None, ax=None, xlim=None, label=None, \
                           show=True):
        '''
        Method to plot a histogram of the integrated pulse areas. Will 
        check to make sure pulses have actually been integrated, plotting
        a spectrum if all conditions are met

        INPUTS

            spectra_type - str, either 'int' or 'ext', refering to spectra 
                of the pulse integrals, or the pulse extrema

            nbin - int, number of bins in the output histogram

            range - 2-tuple of floats, upper and lower limit on pulse area 

            show - boolean, whether to show the plots or just return the
                fig and axes objects

        OUTPUTS

            fig - matplotlib.pyplot.fig, handle for the figure object.
                Useful to overlay stuff on default plots

            ax - matplotlib.pyplot.axes, handle for the axes object 
                containing all the plotted things

        '''

        ### Peform some simple cuts on the data.
        ###    FOR NOW, the only cut peformed is just when the baseline
        ###    is bad, which is determined via a very simple method
        if exclude_bad_baseline:
            good_baseline_inds = np.logical_not(self.bad_baselines)
        else:
            good_baseline_inds = np.ones(self.n_waveform, dtype='bool')

        ### Select the values, either integrals or extrema
        if spectra_type == 'int':
            good_vals = self.waveform_integrals[good_baseline_inds]
        elif spectra_type == 'ext':
            good_vals = self.waveform_extrema[good_baseline_inds]
        else:
            raise ValueError("Argument 'spectra_type' not valid")

        ### If a specific color wasn't chosen, assign the pyplot C0. 
        ### Kind of assuming only one spectra will be plotted at a time
        ### in the naive fashion, and when the user is intentionally 
        ### plotting multiple spectra, they will appropriately use
        ### the color argument
        if color is None:
            color = 'C0'

        ### Generate the figure and axes instances, as well as the value
        ### of the ylimits, since the automatically chosen limits need
        ### to respect ALL the spectra included in the integral
        if fig is None or ax is None:
            fig, ax = plt.subplots(1,1)
            ylim = None
        else:
            ylim = ax.get_ylim()

        if log_scale:
            ax.set_yscale('log')
            if ylim[0] <= 0:
                ylim = (0.5, ylim[1])

        if filled:
            ### If a filled histogram is desired, use the basic pyplot hist
            ### method to obtain the bin values, as well as do the plotting
            bin_vals, bin_edges, _ = \
                ax.hist(good_vals, bins=nbin, range=hist_range, \
                        color=color, label=label)
            ax.set_ylabel('Counts [#]')

        else:
            ### If a non-filled histogram is desired, we probably want to 
            ### normalize to count rate, rather than raw counts, in order
            ### to make comparisons, assuming that's the point of using
            ### a non-filled histogram. Thus we need to do our own plotting

            ### Generate the histogram
            bin_vals, bin_edges = \
                np.histogram(good_vals, bins=nbin, range=hist_range)

            ### Build the array of x values to plot by duplicating the 
            ### bin-edge array, so that vertical sections can be drawn to 
            ### mimic the edge of a histogram bin. Flatten via the FORTRAN 
            ### method to get the right ordering
            plot_x = np.stack((bin_edges, bin_edges)).flatten('F')

            ### Build the corresponding array of y values to plot, by
            ### duplicating the bin values, and then concatenating a half
            ### count to either end, to bring the spectrum back down
            double_bin_vals = np.stack((bin_vals, bin_vals)).flatten('F')
            plot_y = np.concatenate(([0.5], double_bin_vals, [0.5]))

            ### Construct a normalization factor given by the full length
            ### of the exposure associated to this file
            norm_fac = 1.0 / np.max(self.utimes)

            ### Plot the damn thing
            ax.plot(plot_x, plot_y*norm_fac, lw=2, color=color, label=label)
            ax.set_ylabel('Count rate [#/s]')


        ### Setup the x axis label appropriately for the type of histogram
        if spectra_type == 'int':
            ax.set_xlabel('Pulse area [arb.]')
        elif spectra_type == 'ext':
            ax.set_xlabel('Pulse height [arb.]')

        ### Adjust the y limits intelligently, by seeing if the new data 
        ### requires an adjustment (e.g. if it's out of range or something)
        if filled:
            if ylim is None:
                ax.set_ylim(0.5, 1.5*np.max(bin_vals))
            else:
                new_max = 1.5*np.max(bin_vals)
                if ylim[1] > new_max:
                    new_max = ylim[1]
                ax.set_ylim(0.5, new_max)
        else:
            if ylim is None:
                ax.set_ylim(0.5*norm_fac, 1.5*np.max(bin_vals)*norm_fac)
            else:
                new_max = 1.5*np.max(bin_vals)*norm_fac
                new_min = 0.5*norm_fac
                if ylim[1] > new_max:
                    new_max = ylim[1]
                if ylim[0] < new_min:
                    new_min = ylim[0]
                ax.set_ylim(new_min, new_max)

        ### If a desired xlim was provided set the limit to that range. 
        ### Otherwise try the hist_range argument, and then finally the min 
        ### and max of the actual spectra being plotted
        if xlim is not None:
            ax.set_xlim(*xlim)
        elif hist_range is not None:
            ax.set_xlim(*hist_range)
        else:
            ax.set_xlim(np.min(bin_edges), np.max(bin_edges))

        fig.tight_layout()
        if label:
            ax.legend(fontsize=10, loc='upper right')

        ### Show the figure, adding the legend when the show command is 
        ### given, so that all traces are included in the label
        if show:
            plt.show()

        return fig, ax









    def save(self, opt_ext=None, verbose=False):
        '''
        Method to save the class to a pre-defined location for later reloading
        so all of the waveforms don't need to be reprocessed. This is useful
        as the pulse maxima finding algorithm can take a while if the waveforms
        are filtered prior to fitting.

        INPUTS

            opt_ext - str, optional extension to add to the default filename
                if some special analysis was performed

        '''

        ### Make sure the extension has an underscore for intelligible names.
        ### If it doesn't have an underscore, add one
        if opt_ext is not None:
            if opt_ext[0] != '_':
                opt_ext = '_' + opt_ext
        else:
            opt_ext = ''

        ### Extract the parent file name. Probably good to use the os.path
        ### library instead of this derpy solution
        parts = self.fname.split('.')
        savename = parts[0] + opt_ext + '.wfm'

        if verbose:
            print('     ----------------------------------------')
            print('Saving WaveformContainer object... ', end=' ')

        pickle.dump(self, open(savename, 'wb'))

        if verbose:
            print('Done!')
            print('Saved to: ', savename)
            print('     ----------------------------------------')
            sys.stdout.flush()

        return None








    def load(self, opt_ext=None, verbose=False):
        '''
        Method to save the class to a pre-defined location for later reloading
        so all of the waveforms don't need to be reprocessed. This is useful
        as the pulse maxima finding algorithm can take a while if the waveforms
        are filtered prior to fitting.

        INPUTS

            opt_ext - str, optional extension to add to the default filename
                if some special analysis was performed

        '''



        ### Make sure the extension has an underscore for intelligible names.
        ### If it doesn't have an underscore, add one
        if opt_ext is not None:
            if opt_ext[0] != '_':
                opt_ext = '_' + opt_ext
        else:
            opt_ext = ''

        ### Extract the parent file name. Probably good to use the os.path
        ### library instead of this derpy solution
        parts = self.fname.split('.')
        loadname = parts[0] + opt_ext + '.wfm'

        old_class = pickle.load( open(loadname, 'rb') )
        self.__dict__.update(old_class.__dict__)

        if verbose:
            print('     ----------------------------------------')
            print('Loaded WaveformContainer from: ')
            print('          ', loadname)
            print('     ----------------------------------------')
            sys.stdout.flush()

        return None















