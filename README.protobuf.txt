test_compress uses Google's Protocol Buffers (protobuf) tools for
reading and writing configuration data efficiently.

There are three parts to protobuf:
  protoc - the executable that read your .proto file and writes code to
           implement it as objects
  include files - needed when you compile the code that protoc wrote
  libraries - also needed to compile

On Windows

  protoc is available as a precompiled download:
    https://developers.google.com/protocol-buffers/docs/downloads
    See "Protocol Compiler x.x.x binary for windows"

  The include files and libraries must be built from the source.
  Download the full source (same page as the protocol compiler download)
  and read "vsprojects/readme.txt".

On Unix

  # To build it, you'll need autoconf, automake, and aclocal installed.

  ./configure
  make
  make check
  sudo make install

  This will install to /usr/local/lib, /usr/local/include, and
  /usr/local/bin. To install somewhere else:

  make install prefix=/my/other/directory 
