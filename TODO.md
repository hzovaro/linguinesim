############################################################################################
#
#	TO DO:
#	- make a wiki
#	- histogram of shifting-and-stacking errors
#	- window down on star: get some really bright magnitude images so that we can see Airy rings
#	- write a plan: what will go in the paper? 
#	- TODAY: get some results: plots of LI success rate vs. magnitude for a centred star in all 3 bands
#	- how to calculate the SNR in the final image?
#	- shifting-and-stacking xcorr: does it make sense to xcorr with the truth image? because IRL this image is not actually available -- better to xcorr with the first image in the array instead perhaps?
#	- implement more exotic LI techniques 
#	- Check the cryostat temperature: confirm with Rob which dark current value to use. dependence on temperature?
# GALFIT SIMULATIONS
#	- better way of calculating b_n
#	- weird pixelation effect: occurs during call to image.rotate()
#	- modify sersic2D function to take as input z, R_trunc, detector width & detector plate scale instead
#
############################################################################################