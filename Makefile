CXXFLAGS =	-O0 -g3 -Wall -fmessage-length=0 -L -lntl -std=c++0x 

OBJS = md5.o GenericAES.o NTLUtils.o MixingBijections.o WBAES.o WBAESGenerator.o LinearAffineEq.o 

LIBS = -lntl -L/opt/local/lib/ -lboost_iostreams -lboost_serialization -lboost_program_options

TARGET =	MGR_NTL
TARGET01 = test_AEStblSpeed

$(TARGET):  $(OBJS) MGR_NTL.o
	$(CXX) -o $(TARGET) MGR_NTL.o $(OBJS) $(LIBS)
	
$(TARGET01):	$(OBJS) test_AEStblSpeed.o
	$(CXX) -o $(TARGET01) test_AEStblSpeed.o $(OBJS) $(LIBS)

all:	$(TARGET) $(TARGET01)

clean:
	rm -f $(OBJS) $(TARGET) $(TARGET01) MGR_NTL.o test_AEStblSpeed.o
