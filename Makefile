OUTDIR = out

build: joy_libert rinocchio

test: joy_libert
	./$(OUTDIR)/jltest

joy_libert: | out
	g++ -O2 -std=c++11 src/joy_libert_test.cpp src/joy_libert.cpp -o $(OUTDIR)/jltest -pthread -lntl -lgmp -lm

rinocchio: | out
	g++ -O2 -std=c++11 src/rinocchio.cpp src/joy_libert.cpp -o $(OUTDIR)/jltest -pthread -lntl -lgmp -lm

out:
	mkdir -p $(OUTDIR)