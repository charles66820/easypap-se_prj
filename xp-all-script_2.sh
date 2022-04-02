XP_FILES=(ssand-xp-all
ssand-xp-TiledVsLazy
ssand-xp-OmpTiledVsOmpLazy

asand-xp-all
asand-xp-TiledVsLazy
asand-xp-OmpTiledVsOmpLazy)

rm -f xp/*.csv xp_pdf/*.pdf

plots/run-xp-ssand-all.py
plots/run-xp-asand-all.py

for XP_FILE in ${XP_FILES[@]} ; do
  plots/easyplot.py --col schedule --row tile_size --delete iterations -if xp/$XP_FILE.csv -of xp_pdf/$XP_FILE.pdf
done



