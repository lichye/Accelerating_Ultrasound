#compile to run on the host
host: beamform.c solution_check.c
	icc -o beamform beamform.c
	icc -o solution_check solution_check.c
#comiple to run on the Xeon Phi
mic: beamform.c solution_check.c
	icc -mmic -lpthread -o beamform.mic beamform.c -lm
	icc -mmic -lpthread -o solution_check.mic -mmic solution_check.c -lm
clean:
	rm -f *.mic *.o
check: test.sh
	./test.sh