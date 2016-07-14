############################################################################################
#
#	TO DO:
#	- write a plan: what will go in the paper? 
#	- TODAY: get some results: plots of LI success rate vs. magnitude for a centred star in all 3 bands
#	- checked the calculations with Rob - nothing seems amiss but the numbers nonetheless seem quite high
#	- how to calculate the SNR in the final image?
#	- shifting-and-stacking xcorr: does it make sense to xcorr with the truth image? because IRL this image is not actually available -- better to xcorr with the first image in the array instead perhaps?
#	- Remember our goal: what is the minimum magnitude of stars that we can image and still guarantee that the LI technique(s) work(s)?
#		Start off with a single star in the centre of the frame. Get a feel for what works and what doesn't. Run the routine ~100 times and see the average success rate (success = all frames aligned). Then move onto random fields of stars: randomised magnitudes or fixed magnitudes? Does the position of the stars in the image make a big difference? Then after this, move onto HST images. 
#
#	- detector saturation
#	- implement more exotic LI techniques 
#	- Check the cryostat temperature: confirm with Rob which dark current value to use
#	- better way of calculating b_n
#	- weird pixelation effect: occurs during call to image.rotate()
#	- modify sersic2D function to take as input z, R_trunc, detector width & detector plate scale instead
#	- dark current vs. temperature?
#
############################################################################################