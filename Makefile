CC = CC
CFLAGS = -g -Wall -fopenmp -O3
LDFLAGS = -lhdf5

SRC = ./source
TARGET = ns_axion_emission
OBS = nscool

all: $(TARGET)

$(TARGET): $(SRC)/$(TARGET).cpp $(SRC)/$(OBS).o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(SRC)/$(OBS).o: $(SRC)/$(OBS).cpp
	$(CC) $(CFLAGS) -o $@ -c $<

clean:
	$(RM) $(TARGET) $(SRC)/*.o
