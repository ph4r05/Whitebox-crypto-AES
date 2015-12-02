#CXXFLAGS = -O0 -g3 -Wall -fmessage-length=0 -L -lntl -lm -std=c++0x
CXXFLAGS = -O2 -g0 -Wall -fmessage-length=0 -std=c++11 -I/usr/local/include -I/opt/local/include

OBJS = md5.o GenericAES.o NTLUtils.o MixingBijections.o WBAES.o WBAESGenerator.o LinearAffineEq.o BGEAttack.o

LIBS = -lntl -L/usr/local/lib -L/opt/local/lib/ -lboost_iostreams-mt -lboost_serialization-mt -lboost_program_options-mt

TARGET = main
TARGET01 = testing

$(TARGET):  $(OBJS) $(TARGET).o
	$(CXX) -o $(TARGET) $(TARGET).o $(OBJS) $(LIBS)
	
$(TARGET01): $(OBJS) $(TARGET01).o
	$(CXX) -o $(TARGET01) $(TARGET01).o $(OBJS) $(LIBS)

all:	$(TARGET) $(TARGET01)

clean:
	rm -f $(OBJS) $(TARGET) $(TARGET01) $(TARGET).o $(TARGET01).o

