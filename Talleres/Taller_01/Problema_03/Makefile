SJ:
	g++ Sistema_Sol_Jupiter.cpp -o Sistema_Sol_Jupiter.out && ./Sistema_Sol_Jupiter.out
	cd plt && gnuplot SJ.p && cd ../
SJ_r:
	g++ Sistema_Sol_Jupiter_Rotado.cpp -o rotado.out && ./rotado.out
	cd plt && gnuplot SJ_R.p && cd ../
SJT:
	g++ Sistema_Sol_Jupiter_Troyano.cpp -o troyano.out && ./troyano.out
	cd plt && gnuplot SJT.p && cd ../

SJT_p:
	g++ Sistema_Sol_Jupiter_Perturbado_Troyano.cpp -o perturbado.out && ./perturbado.out
	cd plt && gnuplot SJT_p.p && cd ../

SJT_p_periodos:
	g++ Sistema_Sol_Jupiter_Perturbado_Troyano_Periodos.cpp -o perturbado_periodos.out && ./perturbado_periodos.out
	cd plt && gnuplot SJT_p_periodos.p && cd ../
clean:
	rm -f *.out *~ *# *.dat
	clear
