$ENV{'TEXINPUTS'} = '../../third_party/ubcbeamer/ubcbeamer_submodule/tex/latex/ubcbeamer//:' . ($ENV{'TEXINPUTS'} || '');

$pdf_mode = 5;
$pdflatex = 'xelatex -interaction=nonstopmode -file-line-error %O %S';
$dvi_mode = 0;
$postscript_mode = 0;
$clean_ext = 'nav snm';

# Faster edit/compile loop if you do not need SyncTeX.
$synctex = 0;

# Keep generated intermediates; avoids extra churn after cleaning.
$cleanup_includes_cusdep_generated = 0;
