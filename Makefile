CHARMC = charmc -O2 

CXX=$(CHARMC)
CXXFLAGS=$(OPTS)

BINARY=JCN

all: $(BINARY)

JCN_Proj: BSlab.o JCN.o Main.o 
	$(CHARMC) -o $@ Main.o JCN.o BSlab.o  -module CkMulticast -tracemode projections $(OPTS)

$(BINARY): BSlab.o JCN.o Main.o 
	$(CHARMC) -o $@ Main.o JCN.o BSlab.o  -module CkMulticast  $(OPTS)

Main.o: Main.C Main.h Main.decl.h Main.def.h
	$(CHARMC) -o $@ Main.C 

JCN.o: JointContourNet.C JointContourNet.h JCN.decl.h JCN.def.h
	$(CHARMC) -o $@ JointContourNet.C

BSlab.o: BoundSlabs.C BoundSlabs.h BoundSlabs.decl.h BoundSlabs.def.h
	$(CHARMC) -o $@ BoundSlabs.C


%.stamp:% 
	$(CHARMC) $< 

Main.decl.h Main.def.h : Main.ci.stamp
BoundSlabs.decl.h BoundSlabs.def.h: BoundSlabs.ci.stamp
JCN.decl.h JCN.def.h: JointContourNet.ci.stamp 

clean: 
	rm *.o *.decl.h *.def.h 
