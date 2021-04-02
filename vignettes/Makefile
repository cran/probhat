default:
	make probhat.pdf
	make probhatx_kernels.pdf
	make remove_other

probhat.pdf: probhat.tex
	pdflatex probhat.tex

probhatx_kernels.pdf: probhatx_kernels.tex
	pdflatex probhatx_kernels.tex

remove_other:
	rm -f *.tex
	rm -f *.aux
	rm -f *.log
	rm -f temp*
	rm -f Rplots*
