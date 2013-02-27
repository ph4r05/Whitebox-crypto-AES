CXXFLAGS =	-O0 -g -Wall -fmessage-length=0 -L -lntl

OBJS =		MGR_NTL.o

LIBS = -lntl

TARGET =	MGR_NTL

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
