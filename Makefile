CC = CC
CFLAGS = -g -Wall -fopenmp -O3 -std=c++20
LDFLAGS = -lhdf5 -lgsl

SRC = ./source
BUILD_DIR = tmp_build_dir
TARGET = ns_axion_emission
OBS = bremsstrahlung.o utils.o process.o pbf_1s0.o nscool.o ns_axion_emission.o

all: $(TARGET)

$(TARGET): $(addprefix $(BUILD_DIR)/,$(OBS))
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(BUILD_DIR)/%.o: $(SRC)/%.cpp | $(BUILD_DIR)
	$(CC) $(CFLAGS) -o $@ -c $<

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

clean:
	$(RM) -r $(TARGET) $(BUILD_DIR)
