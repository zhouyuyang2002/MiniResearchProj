/**
 * @file CMSketchTest.h
 * @author dromniscience (you@domain.com)
 * @brief Test Count Min Sketch
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once

#include <common/test.h>
#include <sketch/ElasticSketch.h>

#define ES_PARA_PATH "ES.para"
#define ES_TEST_PATH "ES.test"
#define ES_DATA_PATH "ES.data"


#define checkarg(var,name) do{\
  if (!parser.parseConfig())\
    return;\
}while(0)

namespace OmniSketch::Test {

/**
 * @brief Testing class for Count Min Sketch
 *
 */
template <int32_t key_len, typename T, typename hash_t = Hash::AwareHash>
class ElSketchTest : public TestBase<key_len, T> {
  using TestBase<key_len, T>::config_file;

public:
  /**
   * @brief Constructor
   * @details Names from left to right are
   * - show name
   * - config file
   * - path to the node that contains metrics of interest (concatenated with
   * '.')
   */
  ELSketchTest(const std::string_view config_file)
      : TestBase<key_len, T>("Elastic Sketch", config_file, EL_TEST_PATH) {}

  /**
   * @brief Test Elastic Filter
   * @details An overriden method
   */
  void runTest() override;
};

} // namespace OmniSketch::Test

//-----------------------------------------------------------------------------
//
///                        Implementation of templated methods
//
//-----------------------------------------------------------------------------

namespace OmniSketch::Test {

template <int32_t key_len, typename T, typename hash_t>
void CMSketchTest<key_len, T, hash_t>::runTest() {
  /**
   * @brief shorthand for convenience
   *
   */
  using StreamData = Data::StreamData<key_len>;

  /// Part I.
  ///   Parse the config file
  ///
  /// Step i.  First we list the variables to parse, namely:
  ///
  int32_t num_heavy_packet, num_heavy_size;   // sketch config
  int32_t num_light_packet, num_light_size;   // sketch config
  int32_t thre_eject, thre_elephant;          // sketch config

  std::string data_file; // data config
  toml::array arr;       // shortly we will convert it to format
  /// Step ii. Open the config file
  Util::ConfigParser parser(config_file);
  if (!parser.succeed()) {
    return;
  }
  /// Step iii. Set the working node of the parser.
  parser.setWorkingNode(ES_PARA_PATH); // do not forget to to enclose it with braces
  /// Step iv. Parse num_bits and num_hash
  checkarg(num_light_packet, "num_light_packet");
  checkarg(num_heavy_packet, "num_heavy_packet");
  checkarg(num_light_size, "num_light_size");
  checkarg(num_heavy_size, "num_heavy_size");
  checkarg(thre_elephant, "thre_elephant");
  checkarg(thre_eject, "thre_eject");
  /// Step v. Move to the data node
  parser.setWorkingNode(ES_DATA_PATH);
  /// Step vi. Parse data and format
  checkarg(data_file, "data");
  checkarg(arr, "format");
  Data::DataFormat format(arr); // conver from toml::array to Data::DataFormat
  /// [Optional] User-defined rules
  ///
  /// Step vii. Parse Cnt Method.
  std::string method;
  Data::CntMethod cnt_method = Data::InLength;
  checkarg(method, "cnt_method");
  if (!method.compare("InPacket")) {
    cnt_method = Data::InPacket;
  }
  Data::HXMethod hx_method = Data::TopK;
  checkarg(method, "hx_method");
  if (!method.compare("Percentile")) {
    hx_method = Data::Percentile;
  }

  /// Part II.
  ///   Prepare sketch and data
  ///
  /// Step i. Initialize a sketch
  std::unique_ptr<Sketch::SketchBase<key_len, T>> ptr(
    new Sketch::ElasticSketch<key_len, num_heavy_size, num_light_size, T, hash_t>
      (num_heavy_packet, num_light_packet, thre_eject, thre_elephant));
  /// remember that the left ptr must point to the base class in order to call
  /// the methods in it

  /// Step ii. Get ground truth
  ///
  ///       1. read data
  StreamData data(data_file, format); // specify both data file and data format
  if (!data.succeed())
    return;
  Data::GndTruth<key_len, T> gnd_truth;
  gnd_truth.getGroundTruth(data.begin(), data.end(), cnt_method);
  Data::GndTruth<key_len, T> gnd_truth_heavy_hitters;
  gnd_truth_heavy_hitters.getHeavyHitter(gnd_truth, num_heavy_hitter,
                                         hx_method);
  ///       2. [optional] show data info
  fmt::print("DataSet: {:d} records with {:d} keys ({})\n", data.size(),
             gnd_truth.size(), data_file);
  /// Step iii. Insert the samples and then look up all the flows
  ///
  ///        1. update records into the sketch
  this->testUpdate(ptr, data.begin(), data.end(),
                   cnt_method); // metrics of interest are in config file
  ///        2. query for all the flowkeys
  this->testQuery(ptr, gnd_truth); // metrics of interest are in config file
  ///        3. size
  this->testSize(ptr);
  ///        4. heavy
  if (hx_method == Data::TopK) {
    this->testHeavyHitter(
        ptr, gnd_truth_heavy_hitters.min(),
        gnd_truth_heavy_hitters); // metrics of interest are in config file
  } else {
    this->testHeavyHitter(
        ptr, std::floor(gnd_truth.totalValue() * num_heavy_hitter + 1),
        gnd_truth_heavy_hitters); // gnd_truth_heavy_hitter: >, yet HashPipe: >=
  }
  ///        5. show metrics
  this->show();

  return;
}

} // namespace OmniSketch::Test

#undef CM_PARA_PATH
#undef CM_TEST_PATH
#undef CM_DATA_PATH

// Driver instance:
//      AUTHOR: dromniscience
//      CONFIG: sketch_config.toml  # with respect to the `src/` directory
//    TEMPLATE: <13, int32_t, Hash::AwareHash>