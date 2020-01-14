# Fortran compiler
FC ?= gfortran

# original file and revised one
ORG := DC3D.f
REV := DC3D.f90

# download the original file
EXIST = $(shell find . -name $(ORG))
WGET = $(shell wget -O $(ORG) http://www.bosai.go.jp/study/application/dc3d/download/DC3Dfortran.txt)

.PHONY: all test testorg so clean

all: obj

obj:
	$(FC) -c $(REV) -o $(REV).o

org:
	$(if $(EXIST), , $(WGET))
	$(FC) -c $(ORG) -o $(ORG).o

test: obj
	$(FC) test/call.f90  $(REV).o -o xtest
	./xtest

testorg: org
	$(FC) test/call.f  $(ORG).o -o xtest_org
	./xtest_org

so:
	$(FC) $(REV) -o $(REV).so -fPIC -shared
	$(if $(EXIST), , $(WGET))
	$(FC) $(ORG)   -o $(ORG).so   -fPIC -shared

clean:
	-rm -f *.o *.so *.mod x* $(ORG)
