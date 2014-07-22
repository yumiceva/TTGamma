all: ttgamma libttgamma.so

ttgamma: ttgamma3Dict.o EventTree.o BaseScaleFactor.h MuonScaleFactor.h BTagScaleFactor.o BaseSelector.h EvtCleaning.h ElectronSelector.h MuonSelector.h ttgamma3.o runttgamma.C
	g++ -fPIC -o ttgamma `root-config --libs` -I`root-config --incdir` ttgamma3Dict.o EventTree.o BTagScaleFactor.o ttgamma3.o runttgamma.C

libttgamma.so: EventTreeDict.C EventTreeDict.o
	g++ -fPIC -o libttgamma.so `root-config --libs` -I`root-config --incdir` EventTreeDict.o EventTree.o -shared

EventTreeDict.o:
	g++ -fPIC -c -I`root-config --incdir` EventTreeDict.C

EventTreeDict.C:
	rootcint -f EventTreeDict.C -c EventTree.h LinkdefEventTree.h

EventTree.o: EventTree.C EventTree.h
	g++ -fPIC -c -I`root-config --incdir` EventTree.C

BTagScaleFactor.o: BTagScaleFactor.cc BTagScaleFactor.h
	g++ -fPIC -c -I`root-config --incdir` BTagScaleFactor.cc

ttgamma3.o: ttgamma3.C ttgamma3.h
	g++ -fPIC -c -I`root-config --incdir` ttgamma3.C

ttgamma3Dict.o: ttgamma3Dict.C ttgamma3Dict.h
	g++ -fPIC -c -I`root-config --incdir` ttgamma3Dict.C

ttgamma3Dict.C: ttgamma3.h
	rootcint -f ttgamma3Dict.C -c ttgamma3.h Linkdef.h

clean:
	rm -f EventTree.o ttgamma3.o ttgamma3Dict.o ttgamma3Dict.h ttgamma3Dict.C ttgamma libttgamma.so EventTreeDict.o EventTreeDict.h EventTreeDict.C BTagScaleFactor.o
