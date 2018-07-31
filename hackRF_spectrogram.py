from spectrogram_flowgraph import spectrogram_flowgraph
import matplotlib.pyplot as plt
import numpy as np
import scipy
import time
import peakutils

#GLOBAL PARAMS
REMOVE_CLOCK_AND_SIDE_CHANNELS = True
SAMP_RATE = 20.0e6
SIZE_FFT = 2048 # Do not change this value without also changing 'vec_fft_len' to same value in spectrogram_flowgraph.py
N_TO_AVG = int((SAMP_RATE/SIZE_FFT)*1.0)
F_START = 19.0e6 # NOTE: (F_STOP - F_START) should be an integer nultiple of (SAMP_RATE - BW_OVERLAP)
F_STOP = 505.0e6
BW_OVERLAP = 2.0e6 # after making spectra of width SAMP_RATE, this is the total bandwidth to be removed in order to avoid Hack RF roll-off effects from edge of band. BW_OVERLAP/2 is removed from each side of spectra
SAVE_OUTPUT_TO_PDF = False # If you want to save the matplotlib figure being produced

STR_FILE_TEMP = '/Users/josaitis/gnuradio_tutorials/fft_accum_temp/fft_accum_temp_TEST.txt'
STR_FILE_MIN = '/Users/josaitis/gnuradio_tutorials/fft_accum_temp/fft_accum_min.txt'
STR_FILE_MAX = '/Users/josaitis/gnuradio_tutorials/fft_accum_temp/fft_accum_max.txt'
STR_FILE_MEAN = '/Users/josaitis/gnuradio_tutorials/fft_accum_temp/fft_accum_mean.txt'

SIDE_CHAN_TO_DEL=0 #This and the variable below are used when 'REMOVE_CLOCK_AND_SIDE_CHANNELS = True'
CENTER_CHAN_TO_DEL=5


####
ii_first_keep = int( ((BW_OVERLAP/2)/SAMP_RATE)*SIZE_FFT ) #the max and min useable indices of an FFT spectra, after acocunting for removal BW_OVERLAP 
ii_last_keep = int( (1-(BW_OVERLAP/2)/SAMP_RATE)*SIZE_FFT )
print('ii first and last to keep (respectively): '+str(ii_first_keep)+', '+str(ii_last_keep))


###### Functions
def analysis_of_sample():
	file = scipy.fromfile(open(STR_FILE_TEMP), dtype=scipy.float32)
	dat1D = np.array(file)
	dat2D = np.reshape( file,(np.size(file)/SIZE_FFT,SIZE_FFT) ) #order them by every size_fft cols to get rid of edge and clock peaks

	if REMOVE_CLOCK_AND_SIDE_CHANNELS: #set the power at each clock peak and edge (+/- side_chan_to_del channels) to zero
		elims = np.arange(int(SIZE_FFT/2 - CENTER_CHAN_TO_DEL), int(SIZE_FFT/2 +CENTER_CHAN_TO_DEL +1),1) #elimiate the clock peaks and the N channels on each side
		elims = np.concatenate((elims, np.arange(0,int(SIDE_CHAN_TO_DEL)),np.arange(int(SIZE_FFT - SIDE_CHAN_TO_DEL),int(SIZE_FFT))))
		elims = np.unique(elims)

		for ii in elims: 
			dat2D[:,ii]=0

	#link for following line: https://stackoverflow.com/questions/30379311/fast-way-to-take-average-of-every-n-rows-in-a-npy-array
	#note, the power(,2) below is to convert from measured amplitude to measured power. Confirmed in by doing input power tests.
	dat2D_mean = dat2D.transpose().reshape(-1,N_TO_AVG).mean(1).reshape(SIZE_FFT,-1).transpose()# This line avgs every R rows together in a 2D array with C cols
	dat2D_min = dat2D.transpose().reshape(-1,N_TO_AVG).min(1).reshape(SIZE_FFT,-1).transpose()
	dat2D_max = dat2D.transpose().reshape(-1,N_TO_AVG).max(1).reshape(SIZE_FFT,-1).transpose()
	#print('Shape of data after transposes: '+str(np.shape(dat2D_mean)))
	f_handle_mean = open(STR_FILE_MEAN,'ab')
	f_handle_min = open(STR_FILE_MIN,'ab')
	f_handle_max = open(STR_FILE_MAX,'ab')

	np.savetxt(f_handle_mean,np.array(np.ravel(dat2D_mean))[ii_first_keep:ii_last_keep])# Don't save the edges of the spectra, set by BW_OVERLAP, because they can have 'roll-off' effects from the filters in the HackRF. Probably can be calibrated out, this is the lazy way of handling the problem
	np.savetxt(f_handle_min,np.array(np.ravel(dat2D_min))[ii_first_keep:ii_last_keep])
	np.savetxt(f_handle_max,np.array(np.ravel(dat2D_max))[ii_first_keep:ii_last_keep])
	f_handle_mean.close()
	f_handle_min.close()
	f_handle_max.close()
	open(STR_FILE_TEMP, 'w').close() #delete contents of temporary file

def file_flush():
	open(STR_FILE_TEMP, 'w').close()
	open(STR_FILE_MEAN, 'w').close()
	open(STR_FILE_MIN, 'w').close()
	open(STR_FILE_MAX, 'w').close()

def open_and_plot():
	arr_min = np.power(np.loadtxt(STR_FILE_MIN),2)
	arr_max = np.power(np.loadtxt(STR_FILE_MAX),2)
	arr_mean = np.power(np.loadtxt(STR_FILE_MEAN),2)

	x=np.linspace((F_START-(SAMP_RATE- BW_OVERLAP)/2),(F_STOP+(SAMP_RATE - BW_OVERLAP)/2),int(ii_last_keep-ii_first_keep)*n_FFTs)

	plt.semilogy(x/1.e6, np.ravel(arr_mean),'k', label = 'Mean')
	#plt.semilogy(x,np.ravel(arr_min),'b',label='Min')
	#plt.semilogy(x,np.ravel(arr_max),'r',label='Max')
	i_peak = peakutils.indexes(arr_mean, thres=0.3)
	print('Peak found at ('+str(x[i_peak]/1.e6)+' MHz with power of '+ str(arr_mean[i_peak])+')')

	plt.legend()
	if SAVE_OUTPUT_TO_PDF:
		plt.savefig(str('/.../path/....'+str(F_START/1.e6)+'-'+str(F_STOP/1.e6)+'MHz_samp_rate_'+str(SAMP_RATE/1.e6)+'MHz.pdf'))

	plt.xlabel('Frequency (MHz)')
	plt.ylabel('Power (Arbitrary HackRF Units)')

	plt.show()
# Main

n_FFTs = int(((F_STOP-F_START)/(SAMP_RATE- BW_OVERLAP) +1))
freqs = np.arange(F_START,F_STOP + (SAMP_RATE- BW_OVERLAP), SAMP_RATE- BW_OVERLAP)

print('flushing temporary and output files...')
file_flush()
print('done.')

FFT_to_file = spectrogram_flowgraph()
FFT_to_file.set_str_file_out(STR_FILE_TEMP)
FFT_to_file.set_samp_rate(SAMP_RATE)
print('Number of final, averaged, spectra that will be made: '+str(n_FFTs))
print('Number of FFTs to average in order to make one of the above spectra '+str(N_TO_AVG))

for ff in freqs:
	FFT_to_file.set_n_head(N_TO_AVG)
	print('setting center freqeuncy (f_c) to '+str(ff/1.e6)+' MHz...')
	FFT_to_file.set_f_c(ff)
	print('f_c set to: '+str( (FFT_to_file.get_f_c()/1.e6) ) +' MHz')
	FFT_to_file.start()
	FFT_to_file.wait()
	FFT_to_file.stop()
	FFT_to_file.blocks_head_0.reset() #this resets the head block counter, so N_TO_AVG spectra can be output for *each* f_c, not just one f_c
	analysis_of_sample()
	

FFT_to_file.stop()
del(FFT_to_file)

open_and_plot()




