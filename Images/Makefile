# Generate .data files for sizes images

%.data: %.jpg
	java -cp .. WaveletSampleImage $<

all: lenna256.data banff512.data chicago1024.data haleakala2048.data squirrel2048.data hubble4096.data collage8192.data

clean:
	rm -f *.data
