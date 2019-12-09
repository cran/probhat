default:
	make probhat.pdf
	make probhatx_ckernels.pdf
	make remove_other

probhat.pdf: probhat.tex
	pdflatex probhat.tex

probhatx_ckernels.pdf: probhatx_ckernels.tex
	pdflatex probhatx_ckernels.tex

remove_other:
	rm -f *.tex
	rm -f *.aux
	rm -f *.log
	rm -f temp*
	rm -f Rplots*
