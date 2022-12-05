/**
 * @file ElasticSketch.h
 * @author Yuyang Zhou (zhouyuyang2002@stu.pku.edu.com)
 * @brief Elastic Sketch
 * For ElasticSketch itself, I implement enery detail in the paper, including compressing.
 * For tests in Paper, I implement heavy-hitter detection and flow size estimation
 * I don't implement heavy-changer detection && flow-size distribution && test for compressing
 * because they require further change in the bottom of OmniSketch.
 *
 * @copyright Copyright (c) 2022
 */
#pragma once

#include <common/hash.h>
#include <common/sketch.h>

#define lightpacket_t T
#define SKETCH_INF (~((T)0))

namespace OmniSketch::Sketch {

int32_t sketchMin(const int32_t &x, const int32_t &y){
  return x < y ? x : y;
}
int32_t sketchSum(const int32_t &x, const int32_t &y){
  return x + y;
}
template<typename array_t> void allocMem(int dima, int dimb, array_t*** arr){
  array_t* mem = new array_t[dima * dimb];
  *array = new array_t*[dima];
  for (int i = 0; i < dima; i++)
    (*array)[i] = mem + i * dimb;
}
template<typename array_t> void freeMem(array_t*** arr){
  delete[] (*arr)[0];
  delete[] (*arr);
  *arr = NULL;
}

template <int32_t key_len, typename T>
class heavypacket_t{
  public:
    T v_positive;
    T v_light;
    const FlowKey<key_len> &flowkey;
    int16_t ejected;
    int16_t setup;
    heavypacket_t():setup(0){}
    heavypacket_t(T __v_positive, T __v_light
                const FlowKey<key_len> &__flowkey, int16_t __ejected):
                v_positive(__v_positive),v_light(__v_light),
                flowkey(__flowkey), ejected(__ejected), setup(1){}
};
template <int32_t key_len, int32_t h_size, typename T>
class heavybucket_t{
  public:
    heavypacket_t<key_len, T> packet[h_size];
    T v_negative;
    heavypacket_t<key_len, T>* lookup_packet(const FlowKey<key_len> &flowkey, int n_h_packet, int index)const{
      for (int i = 0; i < h_size; i++){
        if (packet[i].flowkey == flowkey)
          return &packet[i];
        else if (packet[i].hashkey % n_h_packet != index)
          return &packet[i];
        else if (packet[i]->set_up == 0)
          return &packet[i];
      }
      int idx = 0;
      T val = ptr[0]->v_positive;
      for (int i = 1; i < h_size; i++)
        if (ptr[i]->v_positive < val){
          val = ptr[i]->v_positive;
          idx = i;
        }
      return &ptr[idx];
    }
}
/**
 * @brief Bloom Filter
 *
 * @tparam key_len  length of flowkey
 * @tparam h_size   size of heavy bucket
 * @tparam l_size   number of hash functions
 * @tparam T        flow size class
 * @tparam hash_t   hashing class
 */
template <int32_t key_len, int32_t h_size, int32_t l_size, typename T, typename hash_t = Hash::AwareHash>
class ElasticSketch : public SketchBase<key_len> {

private:
  int32_t n_h_packet;
  int32_t n_l_packet;
  
  int32_t thre_eject;
  int32_t thre_elephant;
  int32_t thre_n_elephant;
  
  hash_t *hash_h_fns;
  hash_t *hash_l_fns;
  heavybucket_t<key_len, h_size, T> *heavy;
  lightpacket_t **light;
  
  ElasticSketch(const ElasticSketch &) = delete;
  ElasticSketch(ElasticSketch &&) = delete;
  ElasticSketch &operator=(ElasticSketch) = delete;
  
private:
  void realloc_heavy(){
    heavybucket_t<key_len, h_size, T>* new_mem;
    new_mem = new heavybucket_t<key_len, h_size, T>[2 * n_h_packet];
    size_t size_heavy = sizeof(heavybucket_t<key_len, h_size, T>) * n_h_packet;
    memcpy(new_mem, heavy, size_heavy);
    memcpy(new_mem + n_h_packet, heavy, size_heavy);
    delete[] heavy;
    heavy = new_mem;
    n_h_packet *= 2;
    thre_n_elephalt = ((n_h_packet + 2) / 3) * ((h_size + 1) / 2); 
  }
  T queryLightSize(const FlowKey<key_len> &flowkey)const{
    T minimum = SKETCH_INF;
    for (int i = 0; i < l_size; i++){
      int l_index = hash_l_fns[i](flowkey) % n_l_packets;
      minimum = sketchMin(light[i][l_index], minimum);
    }
    return minimum;
  }
  T querySize(const FlowKey<key_len> &flowkey)const{
    int32_t h_index = hash_h_fns[0](flowkey);
    heavypacket_t<key_len, T>* ptr = heavy[h_index]->lookup_packet(flowkey, h_size, h_index);
    if (ptr->setup == 0)
      return 0;
    else{
      int32_t result = 0;
      if (ptr->flowkey == flowkey){
        if (ptr->ejected == 1){
          if (ptr->v_positive + ptr->v_light >= thre_elephant)
            --n_elephant;
          ptr->v_light=queryLightSize(flowkey);
          if (ptr->v_positive + ptr->v_light >= thre_elephant)
            ++n_elephant;
        }
        return ptr->v_positive + ptr->v_light;
      else
        return queryLightSize(flowkey);
    }
  }
}
public:
  /**
   * @brief Construct by specifying # of bits and hash classes
   *
   * @param num_bits        # bit
   * @param num_hash_class  # hash classes
   */
  ElasticSketch(int32_t num_h_packet, int32_t num_l_packet);
  /**
   * @brief Destructor
   *
   */
  ~ElasticSketch();


  /**
   * @brief query the flow size of the key flowkey
   * @details An overriding method
   *
   */
  T query(const FlowKey<key_len> &flowkey) const override;
  /**
   * @brief change the size of flow flowkey, adding the size by value.
   * @details An overriding method
   *
   */
  void update(const FlowKey<key_len> &flowkey, T value) const override;
  /**
   * @brief a quicker version for update
   * @details An overriding method
   *
   */
  void update_quick(const FlowKey<key_len> &flowkey, T value) const override;
  /**
   * @brief look up for a flowkey to see the number of it in the flow.
   * @details An overriding method
   *
   */
  bool lookup(const FlowKey<key_len> &flowkey) const override;
  /**
   * @brief get the heavy hitters, which has flowsize greater than threshold.
   * @details An overriding method
   *
   */
  Data::Estimation<key_len, T> getHeavyHitter(double threshold)const override;
  /**
   * @brief compress the size of the Sketch structure.
   * @details An overriding method
   *
   */
  void compress(int32_t new_num_l_packet) const override;
  /**
   * @brief Size of the sketch
   * @details An overriding method
   */
  size_t size() const override;
  /**
   * @brief Reset the Bloom Filter
   * @details A non-overriding method
   */
  void clear();
};

} // namespace OmniSketch::Sketch

//-----------------------------------------------------------------------------
//
///                        Implementation of templated methods
//
//-----------------------------------------------------------------------------

namespace OmniSketch::Sketch {

template <int32_t key_len, int32_t h_size, int32_t l_size, typename T, typename hash_t = Hash::AwareHash>
ElasticSketch<key_len, T, hash_t>::ElasticSketch(int32_t num_h_packet, int32_t num_l_packet, T __thre_eject, T __thre_elephant){
    :n_h_packet(num_h_packet), n_l_packet(num_l_packet),
    thre_eject(__thre_eject), thre_elephant(__thre_elephant){
  hash_l_fns = new hash_t[l_size];
  hash_h_fns = new hash_t[1];
  allocMem(l_size, n_l_packet, &light);
  heavy = new heavybucket_t<key_len, h_size, T>[n_h_packet];
  thre_n_elephalt = ((n_h_packet + 2) / 3) * ((h_size + 1) / 2); 
}

template <int32_t key_len, int32_t h_size, int32_t l_size, typename T, typename hash_t = Hash::AwareHash>
ElasticSketch<key_len, T, hash_t>::~ElasticSketch(){
  delete[] heavy;
  delete[] hash_l_fns;
  delete[] hash_h_fns;
  freeMem(&light);
}

template <int32_t key_len, int32_t h_size, int32_t l_size, typename T, typename hash_t = Hash::AwareHash>
T ElasticSketch<key_len, T, hash_t>::update(const FlowKey<key_len> &flowkey, T value){
  int32_t h_index = hash_h_fns[0](flowkey) % n_h_packet;
  heavypacket_t<key_len, T>* ptr = heavy[h_index].lookup_packet(flowkey, n_h_packet, h_index);
  if (ptr->setup == 0){
    // not setup
    *ptr = heavypacket_t<key_len, T>(value,0,flowkey,0);
    if (ptr->v_positive + ptr->v_light >= thre_elephant)
      ++n_elephant;
  }
  else if (hash_h_fns[0](ptr->flowkey) % n_h_packet != h_index){
    //fake copies
    //insert is OK, no eject happened
    //reset the v_negative value to  because for heavy resize.
    heavy[h_index].v_negative = 1;
    *ptr = heavypacket_t<key_len, T>(value,queryLightSize(flowkey),flowkey,1);
    if (ptr->v_positive + ptr->v_light >= thre_elephant)
      ++n_elephant;
  }
  else{
    if (ptr->flowkey == flowkey){
      //positive vote
      if (ptr->v_positive + ptr->v_light >= thre_elephant)
        --n_elephant;
      ptr->v_positive += flowkey;
      if (ptr->v_positive + ptr->v_light >= thre_elephant)
        ++n_elephant;
    }
    else{
      heavy[h_index].v_negative += flowkey;
      if (heavy[h_index].v_negative > ptr->v_positive * thre_eject){
        //negative vote: eject current
        int32_t num = ptr->v_positive;
        for (int i=0;i<l_size;i++){
          int l_index = hash_l_fsh[i](ptr->.flowkey)%n_l_packet;
          light[i][l_index] += num;
        }
        if (ptr->v_positive + ptr->v_light >= thre_elephant)
          --n_elephant;
        *ptr = heavypacket_t<keylen, T>(1,queryLightSize(flowkey),flowkey,1);
        heavy[h_index].v_negative = 1;
        if (ptr->v_positive + ptr->v_light >= thre_elephant)
          ++n_elephant;
      }
      else{
        //negative vote: eject insert
        for (int i=0;i<l_size;i++){
          int l_index=hash_l_fsh[i](flowkey)%n_l_packet;
          light[i][l_index] += value;
        }
      }
    }
  }
  if (n_elephant>thre_n_elephant)
    realloc_heavy();
}


template <int32_t key_len, int32_t h_size, int32_t l_size, typename T, typename hash_t = Hash::AwareHash>
void ElasticSketch<key_len, T, hash_t>::update_quick(const FlowKey<key_len> &flowkey, T value){
  int32_t h_index = hash_h_fns[0](flowkey) % n_h_packet;
  heavypacket_t<key_len, T>* ptr = heavy[h_index].lookup_packet(flowkey) n_h_packet;
  if (ptr->setup == 0){
    // not setup
    // ignore elephant check because threshold >> 1
    *ptr = heavypacket_t<key_len, T>(value,0,flowkey,0);
    if (ptr->v_positive + ptr->v_light >= thre_elephant)
      ++n_elephant;
  }
  else if (hash_h_fns[0](ptr->flowkey) % n_h_packet != h_index){
    // fake copies
    // insert is OK, no eject happened
    // ignore Light queries, may result in a fake v_light value
    // but should be fixed in further ejects or queries
    // ignore elephant check because threshold >> 1
    heavy[h_index].v_negative = 1;
    *ptr = heavypacket_t<key_len, T>(value, 0, flowkey, 1);
    if (ptr->v_positive + ptr->v_light >= thre_elephant)
      ++n_elephant;
  }
  else{
    if (ptr->flowkey == flowkey){
      if (ptr->v_positive + ptr->v_light == thre_elephant)
        --n_elephant;
      //positive vote
      ptr->v_positive += value;
      if (ptr->v_positive + ptr->v_light == thre_elephant)
        ++n_elephant;
    }
    else{
      ptr->v_negative += value;
      if (ptr->v_negative > ptr->v_positive * thre_eject){
        //negative vote: eject current
        //ignore light update
        if (ptr->v_positive + ptr->v_light >= thre_elephant)
          --n_elephant;
        //ignore light queries
        heavy[h_index].v_negative = 1;
        *ptr = heavypacket_t<key_len, T>(ptr->v_positive + 1,1,0,flowkey,1);
        // ignore elephant check because threshold >> 1
        if (ptr->v_positive + ptr->v_light >= thre_elephant)
          ++n_elephant;
      }
      else{
        //negative vote: eject insert
        //ignore light update
      }
    }
  }
  //ignore heavy realloc
}

template <int32_t key_len, int32_t h_size, int32_t l_size, typename T, typename hash_t = Hash::AwareHash>
bool ElasticSketch<key_len, T, hash_t>::lookup(const FlowKey<key_len> &flowkey)const{
  return querySize(flowkey) != 0;
}


template <int32_t key_len, int32_t h_size, int32_t l_size, typename T, typename hash_t = Hash::AwareHash>
T ElasticSketch<key_len, T, hash_t>::query(const FlowKey<key_len> &flowkey)const{
  return querySize(flowkey);
}


template <int32_t key_len, int32_t h_size, int32_t l_size, typename T, typename hash_t = Hash::AwareHash>
Data::Estimation<key_len, T> ElasticSketch< key_len, T >::getHeavyHitter(double threshold)const{
  Data::Estimation<key_len, T> result;
  for (int i = 0; i < n_h_packet; i++)
    for (int j = 0; j < h_size; j++)
      if (heavy[i].packet[j].setup == 0)
        break;
      else{
        T value = querySize(heavy[i].packet[j].flowkey);
        if (value >= threshold)
          result[heavy[i].packet[j].flowkey] = value;
      }
  return result;
}

template <int32_t key_len, int32_t h_size, int32_t l_size, typename T, typename hash_t = Hash::AwareHash>
int32_t ElasticSketch<key_len, hash_t>::compress(int32_t new_num_l_packet, bool method){
  if (n_l_packet % new_num_l_packet != 0 || new_num_l_packet != 0)
    return -1; // error
  if (n_l_packet == new_num_l_packet)
    return 0;
  int32_t (*func)(int32_t, int32_t) = (method ? sketchMax : sketchSum);
  lightpacket_t** new_mem;
  allocMem(l_size, new_num_l_packet);
  for (int i = 0; i < l_size; i++){
    lightpacket_t* arr_old = light[i];
    lightpacket_t* arr_new = new_mem[i];
    memcpy(arr_new, arr_old, sizeof(lightpacket_t) * new_num_l_packet);
    int round = n_l_packet / new_num_l_packet;
    for (int j = 1; j < round; j++){
      lightpacket_t* ptr_old = &arr_old[j * new_num_l_packet];
      lightpacket_t* ptr_new = arr_new;
      for (int k = 0; k < new_num_l_packet; k++, ptr_old++, ptr_new++)
        *ptr_new = func(*ptr_new, *ptr_old);
    }
  }
  n_l_packet = new_num_l_packet;
  freeMem(&light)
  light = new_mem;
  return 0;
}

template <int32_t key_len, int32_t h_size, int32_t l_size, typename T, typename hash_t = Hash::AwareHash>
size_t ElasticSketch<key_len, hash_t>::size() const {
  return sizeof(*this)                // Instance
         + n_h_packet * sizeof(heavybucket_t<key_len, size, T>)  // heavy
         + l_hash * sizeof(lightpacket_t*)    // light
		 + l_size * n_l_packet * sizeof(lightpacket_t) //light[0]
}

template <int32_t key_len, int32_t h_size, int32_t l_size, typename T, typename hash_t = Hash::AwareHash>
void ElasticSketch<key_len, hash_t>::clear(){
  std::fill(heavy[0], heavy[0] + n_h_packet * sizeof(heavybucket_t<key_len, h_size, T>), 0);
  std::fill(light[0], light[0] + l_size * n_l_packet * sizeof(lightpacket_t), 0);
}

} // namespace OmniSketch::Sketch

#undef lightpacket_t
#undef SKETCH_INF