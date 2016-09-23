CHARMC = charmc -std=c++0x -O2 

CXX=$(CHARMC)
CXXFLAGS=$(OPTS)

BINARY=JCN

all: $(BINARY)

JCN_Proj: BSlab.o JCN.o Main.o  
	$(CHARMC) -o $@ Main.o JCN.o BSlab.o RectGData.o CompJCN.o MrgJCN.o PtopeGeom.o PyramidTree.o -module CkMulticast -tracemode projections $(OPTS)

$(BINARY): RectGData.o BSlab.o JCN.o Main.o 
	$(CHARMC) -o $@ Main.o JCN.o BSlab.o RectGData.o CompJCN.o MrgJCN.o PtopeGeom.o PyramidTree.o -module CkMulticast $(OPTS)

Main.o: Main.C Main.h Main.decl.h Main.def.h
	$(CHARMC) -o $@ Main.C 

JCN.o: RectGData.o CompJCN.o MrgJCN.o JointContourNet.C JointContourNet.h JCN.decl.h JCN.def.h  
	$(CHARMC) -o $@ JointContourNet.C 

BSlab.o: BoundSlabs.C BoundSlabs.h BoundSlabs.decl.h BoundSlabs.def.h CommonDef.h
	$(CHARMC) -o $@ BoundSlabs.C

RectGData.o:  RectilinearGridData.C RectilinearGridData.h CommonDef.h
	$(CHARMC) -o $@ RectilinearGridData.C 

CompJCN.o: PtopeGeom.o PyramidTree.o CompJointContourNet.C CompJointContourNet.h 
	$(CHARMC) -o $@ CompJointContourNet.C

MrgJCN.o: MergJointContourNet.C MergJointContourNet.h
	$(CHARMC) -o $@ MergJointContourNet.C

PtopeGeom.o: PolytopeGeometry.C PolytopeGeometry.h 
	$(CHARMC) -o $@ PolytopeGeometry.C

PyramidTree.o: PyramidTree.C PyramidTree.h CommonDef.h
	$(CHARMC) -o $@ PyramidTree.C

%.stamp:% 
	$(CHARMC) $< 

Main.decl.h Main.def.h : Main.ci.stamp
BoundSlabs.decl.h BoundSlabs.def.h: BoundSlabs.ci.stamp
JCN.decl.h JCN.def.h: JointContourNet.ci.stamp 

clean: 
	rm *.o *.decl.h *.def.h 
