CC = CC
CFLAGS = -g -Wall -fopenmp -O3
LDFLAGS = -lhdf5

SRC = ./source
TARGET = ns_axion_emission
OBS = pbf_process.o nscool.o

all: $(TARGET)

$(TARGET): $(SRC)/$(TARGET).cpp $(addprefix $(SRC)/,$(OBS))
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(SRC)/%.o: $(SRC)/%.cpp
	$(CC) $(CFLAGS) -o $@ -c $<

clean:
	$(RM) $(TARGET) $(addprefix $(SRC)/,$(OBS))
