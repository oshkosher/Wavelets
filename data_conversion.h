

class DataConversion {

  enum Type {
    INT8, INT32, FLOAT32, FLOAT64
  };

  setFromType(Type t);
  setToType(Type t);

  class Range {
    setAdaptive();  // use the minimum and maximum values found in the data
    set(double min, double max);  // inclusive
  };

  Range fromRange, toRange;

};
