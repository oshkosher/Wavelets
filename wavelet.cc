#include "wavelet.h"

QuantizeAlgorithm quantAlgName2Id(const char *name) {
  if (!strcmp(name, "uniform")) return QUANT_ALG_UNIFORM;
  if (!strcmp(name, "log")) return QUANT_ALG_LOG;
  if (!strcmp(name, "lloyd")) return QUANT_ALG_LLOYD;
  return QUANT_ALG_UNKNOWN;
}
  

const char *quantAlgId2Name(QuantizeAlgorithm id) {
  switch (id) {
  case QUANT_ALG_UNIFORM: return "uniform";
  case QUANT_ALG_LOG: return "log";
  case QUANT_ALG_LLOYD: return "lloyd";
  default: return NULL;
  }
}
  

CompressionParametersBuffer_QuantizationAlgorithm quantAlgId2ProtoId
  (QuantizeAlgorithm id) {

  switch (id) {
  case QUANT_ALG_UNIFORM:
    return CompressionParametersBuffer_QuantizationAlgorithm_UNIFORM;
  case QUANT_ALG_LOG: 
    return CompressionParametersBuffer_QuantizationAlgorithm_LOG;
  case QUANT_ALG_LLOYD:
    return CompressionParametersBuffer_QuantizationAlgorithm_LLOYD;
  default:
    assert(false);
    return CompressionParametersBuffer_QuantizationAlgorithm_LOG;
  }
}
  

QuantizeAlgorithm quantProtoId2AlgId
  (CompressionParametersBuffer_QuantizationAlgorithm protoId) {

  switch (protoId) {
  case CompressionParametersBuffer_QuantizationAlgorithm_UNIFORM:
    return QUANT_ALG_UNIFORM;
  case CompressionParametersBuffer_QuantizationAlgorithm_LOG:
    return QUANT_ALG_LOG;
  case CompressionParametersBuffer_QuantizationAlgorithm_LLOYD:
    return QUANT_ALG_LLOYD;
  default:
    return QUANT_ALG_UNKNOWN;
  }
}

// convert between my wavelet algorithm enum and the protobuf enum
CompressionParametersBuffer_WaveletAlgorithm waveletAlgToProtoId
  (WaveletAlgorithm id) {
  switch (id) {
  case WAVELET_HAAR:
    return CompressionParametersBuffer_WaveletAlgorithm_HAAR;
  case WAVELET_CDF97:
    return CompressionParametersBuffer_WaveletAlgorithm_CDF97;
  case WAVELET_CDF53:
    return CompressionParametersBuffer_WaveletAlgorithm_CDF53;
  default:
    assert(false);
    return CompressionParametersBuffer_WaveletAlgorithm_CDF97;
  }
}

WaveletAlgorithm protoIdToWaveletAlg
  (CompressionParametersBuffer_WaveletAlgorithm protoId) {
  switch (protoId) {
  case CompressionParametersBuffer_WaveletAlgorithm_HAAR:
    return WAVELET_HAAR;
  case CompressionParametersBuffer_WaveletAlgorithm_CDF97:
    return WAVELET_CDF97;
  case CompressionParametersBuffer_WaveletAlgorithm_CDF53:
    return WAVELET_CDF53;
  default:
    assert(false);
    return WAVELET_UNKNOWN;
  }
}


WaveletAlgorithm waveletAlgNameToId(const char *name) {
  if (!strcmp(name, waveletAlgToName(WAVELET_HAAR))) return WAVELET_HAAR;
  if (!strcmp(name, waveletAlgToName(WAVELET_CDF97))) return WAVELET_CDF97;
  if (!strcmp(name, waveletAlgToName(WAVELET_CDF53))) return WAVELET_CDF53;
  return WAVELET_UNKNOWN;
}


const char *waveletAlgToName(WaveletAlgorithm id) {
  switch (id) {
  case WAVELET_HAAR: return "haar";
  case WAVELET_CDF97: return "cdf97";
  case WAVELET_CDF53: return "cdf53";
  default: return "unknown";
  }
}

// convert between my wavelet data type enum and the protobuf enum
CubeletBuffer_DataType datatypeToProtoId(WaveletDataType id) {
  switch (id) {
  case WAVELET_DATA_UINT8:
    return CubeletBuffer_DataType_UINT8;
  case WAVELET_DATA_INT32:
    return CubeletBuffer_DataType_INT32;
  case WAVELET_DATA_FLOAT32:
    return CubeletBuffer_DataType_FLOAT32;
  default:
    assert(false);
    return CubeletBuffer_DataType_UINT8;
  }
}

WaveletDataType protoIdToWaveletDatatype(CubeletBuffer_DataType id) {
  switch (id) {
  case CubeletBuffer_DataType_UINT8:
    return WAVELET_DATA_UINT8;
  case CubeletBuffer_DataType_INT32:
    return WAVELET_DATA_INT32;
  case CubeletBuffer_DataType_FLOAT32:
    return WAVELET_DATA_FLOAT32;
  default:
    assert(false);
    return WAVELET_DATA_UNKNOWN;
  }
}

const char *waveletDataTypeName(WaveletDataType id) {
  switch (id) {
  case WAVELET_DATA_UINT8: return "UINT8";
  case WAVELET_DATA_INT32: return "INT32";
  case WAVELET_DATA_FLOAT32: return "FLOAT32";
  default: return "UNKNOWN";
  }
}

int waveletDataTypeSize(WaveletDataType id) {
  switch (id) {
  case WAVELET_DATA_UINT8: return 1;
  case WAVELET_DATA_INT32: 
  case WAVELET_DATA_FLOAT32: return 4;
  default: return -1;
  }
}

template<> void CubeNum<float>::setType() {
  datatype = WAVELET_DATA_FLOAT32;
}

template<> const char *CubeNum<float>::printfFormat() {
  return "%8.4f ";
}

/*
template<> void CubeNum<float>::print() {
  // for (int z=0; z<depth(); z++) {
  for (int z=0; z<1; z++) {
    printf("z=%d\n", z);
    for (int y=0; y<height(); y++) {
      for (int x=0; x<width(); x++) {
        printf("%8.4f ", get(x,y,z));
      }
      putchar('\n');
    }
    putchar('\n');
  }
}
*/

template<> void CubeNum<int>::setType() {
  datatype = WAVELET_DATA_INT32;
}

template<> const char *CubeNum<int>::printfFormat() {
  return "%5d ";
}

/*
template<> void CubeNum<int>::print() {
  // for (int z=0; z<depth(); z++) {
  for (int z=0; z<1; z++) {
    printf("z=%d\n", z);
    for (int y=0; y<height(); y++) {
      for (int x=0; x<width(); x++) {
        printf("%5d ", get(x,y,z));
      }
      putchar('\n');
    }
    putchar('\n');
  }
}
*/

template<> void CubeNum<unsigned char>::setType() {
  datatype = WAVELET_DATA_UINT8;
}

template<> const char *CubeNum<unsigned char>::printfFormat() {
  return "%3d ";
}

/*
template<> void CubeNum<unsigned char>::print() {
  // for (int z=0; z<depth(); z++) {
  for (int z=0; z<1; z++) {
    printf("z=%d\n", z);
    for (int y=0; y<height(); y++) {
      for (int x=0; x<width(); x++) {
        printf("%3d ", get(x,y,z));
      }
      putchar('\n');
    }
    putchar('\n');
  }
}
*/


void Cube::copyFromCubeletBuffer(const CubeletBuffer *buf) {
  init();

  size = int3(buf->width(), buf->height(), buf->depth());
  // cubelet buffer doesn't encode an inset, so leave it initialized to 0
  totalSize = size;
  parentOffset = int3(buf->x_offset(), buf->y_offset(), buf->z_offset());

  dataFileOffset = buf->data_file_offset();
  datatype = protoIdToWaveletDatatype(buf->data_type());
  if (buf->compression_algorithm() !=
      CubeletBuffer_CompressionAlgorithm_WAVELET) {
    isWaveletCompressed = false;
  } else {
    isWaveletCompressed = true;
    param.compressedSize = buf->byte_count();

    const CompressionParametersBuffer &p = buf->compression_parameters();
    param.originalSize = int3
      (p.original_width(),
       p.original_height(),
       p.original_depth());
    param.originalDatatype = protoIdToWaveletDatatype(p.original_data_type());
    param.transformSteps = int3
      (p.wavelet_transform_steps_x(),
       p.wavelet_transform_steps_y(),
       p.wavelet_transform_steps_z());
    param.waveletAlg = protoIdToWaveletAlg(p.wavelet_algorithm());
    param.isWaveletTransposeStandard = p.standard_transpose();
    param.thresholdFraction = p.threshold_fraction();
    param.thresholdValue = p.threshold_value();
    param.binCount = p.quant_bin_count();
    param.quantAlg = quantProtoId2AlgId(p.quantization_algorithm());

    if (param.quantAlg == QUANT_ALG_UNIFORM ||
        param.quantAlg == QUANT_ALG_LOG) {
      param.maxValue = p.quant_max_value();
    } else {
      param.binBoundaries.resize(p.quant_bin_boundaries_size());
      for (int i=0; i < p.quant_bin_boundaries_size(); i++)
        param.binBoundaries[i] = p.quant_bin_boundaries(i);

      param.binValues.resize(p.quant_bin_values_size());
      for (int i=0; i < p.quant_bin_values_size(); i++)
        param.binValues[i] = p.quant_bin_values(i);
    }

    param.huffDecode.resize(p.huffman_decode_table_size());
    for (int i=0; i < p.huffman_decode_table_size(); i++)
      param.huffDecode[i] = p.huffman_decode_table(i);
    
  }
}

void Cube::copyToCubeletBuffer(CubeletBuffer *buf) const {
  buf->set_width(width());
  buf->set_height(height());
  buf->set_depth(depth());

  buf->set_x_offset(parentOffset.x);
  buf->set_y_offset(parentOffset.y);
  buf->set_z_offset(parentOffset.z);

  if (dataFileOffset > 0) buf->set_data_file_offset(dataFileOffset);
  buf->set_data_type(datatypeToProtoId(datatype));

  if (!isWaveletCompressed) {
    buf->set_compression_algorithm(CubeletBuffer_CompressionAlgorithm_NONE);
  } else {
    buf->set_compression_algorithm(CubeletBuffer_CompressionAlgorithm_WAVELET);
    buf->set_byte_count(param.compressedSize);

    CompressionParametersBuffer *p = buf->mutable_compression_parameters();
    p->set_original_width(param.originalSize.x);
    p->set_original_height(param.originalSize.y);
    p->set_original_depth(param.originalSize.z);

    if (param.originalDatatype != WAVELET_DATA_UNKNOWN)
      p->set_original_data_type(datatypeToProtoId(param.originalDatatype));
    
    p->set_wavelet_transform_steps_x(param.transformSteps.x);
    p->set_wavelet_transform_steps_y(param.transformSteps.y);
    p->set_wavelet_transform_steps_z(param.transformSteps.z);

    p->set_wavelet_algorithm(waveletAlgToProtoId(param.waveletAlg));
    p->set_standard_transpose(param.isWaveletTransposeStandard);
    p->set_threshold_fraction(param.thresholdFraction);
    p->set_threshold_value(param.thresholdValue);
    p->set_quant_bin_count(param.binCount);
    p->set_quantization_algorithm(quantAlgId2ProtoId(param.quantAlg));

    if (param.quantAlg == QUANT_ALG_UNIFORM ||
        param.quantAlg == QUANT_ALG_LOG) {
      p->set_quant_max_value(param.maxValue);
    } else {
      for (size_t i=0; i < param.binBoundaries.size(); i++)
        p->add_quant_bin_boundaries(param.binBoundaries[i]);
      for (size_t i=0; i < param.binValues.size(); i++)
        p->add_quant_bin_values(param.binValues[i]);
    }

    for (size_t i=0; i < param.huffDecode.size(); i++)
      p->add_huffman_decode_table(param.huffDecode[i]);

  }
}

  

