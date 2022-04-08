rm -f xp/*.csv xp_pdf/*.pdf

plots/run-xp-ssand-3.py

plots/easyplot.py --col tiling --delete iterations -if xp/ssand-xp-avx.csv -of xp_pdf/ssand-xp-avx.pdf
#plots/easyplot.py -if xp/ssand-xp-ocl-perf.csv --plottype heatmap -heatx tilew -heaty tileh -of xp_pdf/ssand-xp-ocl-perf.pdf
# plots/easyplot.py --col tiling --row tile_size --delete iterations -if xp/sand-xp-ocl.csv -of xp_pdf/sand-xp-ocl.pdf



