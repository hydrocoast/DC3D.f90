# Fortran compiler
FC ?= gfortran

# original file and revised one
ORG := DC3D.f
REV := DC3D.f90

# download the original file
EXIST = $(shell find . -name $(ORG))
WGET = $(shell wget -O $(ORG) https://www.bosai.go.jp/information/pdf/DC3Dfortran.txt)

.PHONY: all clean

all: obj

obj:
	$(FC) -c $(REV) -o $(REV).o

getorg:
	$(if $(EXIST), , $(WGET))

org: getorg
	$(FC) -c $(ORG) -o $(ORG).o

test: obj
	$(FC) test/call.f90  $(REV).o -o xtest
	./xtest

testorg: org
	$(FC) test/call.f  $(ORG).o -o xtest_org
	./xtest_org

so: getorg
	$(FC) $(REV) -o $(REV).so -fPIC -shared
	$(FC) $(ORG)   -o $(ORG).so   -fPIC -shared

clean:
	-rm -f *.o *.so *.mod x* $(ORG)
