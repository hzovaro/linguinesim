############################################################################################

Organisational stuff:
- Make a Wiki - discuss preferred format with Francis. Preferably somewhere with LaTeX formatting
- LaTeX template for paper w/ section headings 
- Make short-term and long-term plans
- Astro computing course?

############################################################################################

Coding stuff:
- Top-level work/data flow
- Decide whether to use MATLAB or Python for phase screen simulations: look at Mike's code (1 hr)
- Research more lucky imaging methods
- histogram of shifting-and-stacking errors
- get some results: plots of LI success rate vs. magnitude for a centred star in all 3 bands
- window down on star: get some really bright magnitude images so that we can see Airy rings

############################################################################################

Galfit simulations:
- better way of calculating b_n
- weird pixelation effect: occurs during call to image.rotate()
- modify sersic2D function to take as input z, R_trunc, detector width & detector plate scale instead

############################################################################################
 
Things to ask:
- Check the cryostat temperature: confirm with Rob which dark current value to use. dependence on temperature?
- shifting-and-stacking xcorr: does it make sense to xcorr with the truth image? because IRL this image is not actually available -- better to xcorr with the first image in the array instead perhaps?
- how to calculate the SNR in the final image?


