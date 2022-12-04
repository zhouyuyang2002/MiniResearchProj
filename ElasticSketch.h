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

#define lightpacket_t int32_t
#define SKETCH_INF (0x3fffffff)

namespace OmniSketch::Sketch {

int32_t sketchMin(const int32_t &x, const int32_t &y){
  return x < y ? x : y;
}
int32_t sketchSum(const int32_t &x, const int32_t &y){
  return x + y;
}
template <int32_t key_len>
class heavypacket_t{
  public:
    int32_t v_positive;
    int32_t v_negative;
    int32_t v_light;
    const FlowKey<key_len> &flowkey;
    int16_t ejected;
    int16_t setup;
    heavypacket_t():setup(0){}
    heavypacket_t(int32_t __v_positive, int32_t __v_negative, int32_t __v_light
                const FlowKey<key_len> &__flowkey, int32_t __ejected):
                v_positive(__v_positive),v_negative(__v_negative),v_light(__v_light),
                flowkey(__flowkey), ejected(__ejected), setup(1){}
};
/**
 * @brief Bloom Filter
 *
 * @tparam key_len  length of flowkey
 * @tparam hash_t   hashing class
 */
template <int32_t key_len, typename T, typename hash_t = Hash::AwareHash>
class ElasticSketch : public SketchBase<key_len> {

private:
  
  int32_t n_h_packet;
  int32_t n_h_size;
  int32_t n_l_packet;
  int32_t n_l_hash;
  
  int32_t thre_eject;
  int32_t thre_elephant;
  int32_t thre_n_elephant;
  
  hash_t *hash_h_fns;
  hash_t *hash_l_fns;
  heavypacket_t<key_len> **heavy;
  T **heavy_val;
  lightpacket_t **light;
  T **light_val;
  
  ElasticSketch(const ElasticSketch &) = delete;
  ElasticSketch(ElasticSketch &&) = delete;
  ElasticSketch &operator=(ElasticSketch) = delete;
  
private:
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
  
  int32_t queryLightSize(const FlowKey<key_len> &flowkey)const{
    int minimum = SKETCH_INF;
    for (int i = 0; i < n_l_hash; i++){
      int l_index = hash_l_fns[i](flowkey)%n_l_packets;
      minimum = sketchMin(light[i][l_index], minimum);
    }
    return minimum;
  }
  
  heavypacket_t* lookup_packet(const FlowKey<key_len> &flowkey, heavypacket_t* ptr)const{
    for (int i = 0; i < n_h_size; i++){
      if (ptr[i]->flowkey == flowkey)
        return &ptr[i];
      else if (ptr[i]->set_up == 0)
        return &ptr[i];
    }
    int idx = 0;
    int val = ptr[0]->v_positive;
    for (int i = 1; i < n_h_size; i++)
      if (ptr[i]->v_positive < val){
        val = ptr[i]->v_positive;
        idx = i;
      }
    return &ptr[idx];
  }
  void realloc_heavy(){
    heavypacket_t** new_mem;
    allocMem(n_h_packet * 2, n_h_size, &new_mem);
    memcpy(new_mem, heavy, sizeof(heavypacket_t) * n_h_packet * n_h_size);
    memcpy(new_mem + n_h_packet, heavy, sizeof(heavypacket_t) * n_h_packet * n_h_size);
    freeMem(&heavy);
    heavy = new_mem;
    n_h_packet *= 2;
    thre_n_elephalt = ((n_h_packet + 2) / 3) * ((n_h_size + 1) / 2); 
  }
  int32_t querySize(const FlowKey<key_len> &flowkey)const{
    int32_t h_index = hash_h_fns[0](flowkey);
    heavypacket_t* ptr = lookup_packet(flowkey, heavy[h_index]);
    if (ptr->setup == 0)
      return 0;
    else{
      int32_t result = 0;
      if (ptr->flowkey == flowkey){
        if (heavy[h_index].ejected == 1){
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
  ElasticSketch(int32_t num_h_packet, int32_t num_h_dup, int32_t num_l_packet, int32_t num_l_hash);
  /**
   * @brief Destructor
   *
   */
  ~ElasticSketch();

  /**
   * @brief Insert a flowkey into the bloom filter
   * @details An overriding method
   *
   */
  void insert(const FlowKey<key_len> &flowkey) override;
  /**
   * @brief Realloc the memory of array heavy, double its size
   * @details An overriding method
   *
   */
  void realloc_heavy() override;
  /**
   * @brief query a flowkey to see the number of it in the light backet.
   * @details An overriding method
   *
   */
  int32_t queryLightSize(const FlowKey<key_len> &flowkey) override;
  /**
   * @brief query a flowkey to see the number of it in the flow.
   * @details An overriding method
   *
   */
  int querySize(const FlowKey<key_len> &flowkey) const override;
  /**
   * @brief look up for a flowkey to see the number of it in the flow.
   * @details An overriding method
   *
   */
  int lookup(const FlowKey<key_len> &flowkey) const override;
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

template <int32_t key_len, typename T, typename hash_t>
ElasticSketch<key_len, T, hash_t>::ElasticSketch(int32_t num_h_packet, int32_t num_h_size, int32_t num_l_packet, int32_t num_l_hash)
    :n_h_packet(num_h_packet), n_h_size(num_h_size), n_l_packet(num_l_packet), n_l_hash(num_l_hash),
    thre_eject(8), thre_elephant(5000), thre_n_elephant((n_h_packet+1)/2), n_elephant(0){
  hash_l_fns = new hash_t[n_l_hash];
  hash_h_fns = new hash_t[1];
  allocMem(n_h_packet, n_h_size, &heavy);
  allocMem(n_l_hash, n_l_packet, &light);
  thre_n_elephalt = ((n_h_packet + 2) / 3) * ((n_h_size + 1) / 2); 
}

template <int32_t key_len, typename T, typename hash_t>
ElasticSketch<key_len, T, hash_t>::~ElasticSketch(){
  freeMem(&heavy);
  freeMem(&light);
}

template <int32_t key_len, typename T, typename hash_t>
T ElasticSketch<key_len, T, hash_t>::update(const FlowKey<key_len> &flowkey, T value){
  
  int32_t h_index = hash_h_fns[0](flowkey);
  heavypacket_t ptr = lookup_packet(flowkey, heavy[h_index]);
  if (ptr->setup == 0){
    // not setup
    *ptr = heavypacket_t<key_len>(value,0,0,flowkey,0);
    if (ptr->v_positive + ptr->v_light >= thre_elephant)
      ++n_elephant;
  }
  else if (hash_h_fns[0](ptr->flowkey)!=h_index){
    //fake copies
    //insert is OK, no eject happened
    *ptr = heavypacket_t(value,1,queryLightSize(flowkey),flowkey,1);
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
      ptr->v_negative += flowkey;
      if (ptr->v_negative > ptr->v_positive * thre_eject){
        //negative vote: eject current
        int32_t num = ptr->v_positive;
        for (int i=0;i<n_l_hash;i++){
          int l_index = hash_l_fsh[i](ptr->.flowkey)%n_l_packet;
          light[i][l_index] += num;
        }
        if (ptr->v_positive + ptr->v_light >= thre_elephant)
          --n_elephant;
        *ptr = heavypacket_t(1,1,queryLightSize(flowkey),flowkey,1);
        if (ptr->v_positive + ptr->v_light >= thre_elephant)
          ++n_elephant;
      }
      else{
        //negative vote: eject insert
        for (int i=0;i<n_l_hash;i++){
          int l_index=hash_l_fsh[i](flowkey)%n_l_packet;
          light[i][l_index] += value;
        }
      }
    }
  }
  if (n_elephant>thre_n_elephant)
    realloc_heavy();
}


template <int32_t key_len, typename T, typename hash_t>
void ElasticSketch<key_len, T, hash_t>::update_quick(const FlowKey<key_len> &flowkey, T value){
  int32_t h_index = hash_h_fns[0](flowkey);
  heavypacket_t ptr = lookup_packet(flowkey, heavy[h_index]);
  if (ptr->setup == 0){
    // not setup
    // ignore elephant check because threshold >> 1
    *ptr = heavypacket_t<key_len>(value,0,0,flowkey,0);
    if (ptr->v_positive + ptr->v_light >= thre_elephant)
      ++n_elephant;
  }
  else if (hash_h_fns[0](ptr->flowkey) != h_index){
    // fake copies
    // insert is OK, no eject happened
    // ignore Light queries, may result in a fake v_light value
    // but should be fixed in further ejects or queries
    // ignore elephant check because threshold >> 1
    *ptr = heavypacket_t(value,1,0,flowkey,1);
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
        *ptr = heavypacket_t(ptr->v_positive + 1,1,0,flowkey,1);
        // ignore elephant check because threshold >> 1
      }
      else{
        //negative vote: eject insert
        //ignore light update
      }
    }
  }
  //ignore heavy realloc
}

template <int32_t key_len, typename T, typename hash_t>
bool ElasticSketch<key_len, T, hash_t>::lookup(const FlowKey<key_len> &flowkey)const{
  return querySize(flowkey) != 0;
}


template <int32_t key_len, typename T, typename hash_t>
OmniSketch::Data::Estimation<key_len, T> MySketch< key_len, T >::getHeavyHitter(double threshold)const{
  Data::Estimation<key_len, T> result;
  for (int i = 0; i < n_h_packet; i++)
    for (int j = 0; j < n_h_size; j++)
      if (heavy[i][j].setup == 0)
        break;
      else{
        T value = query(heavy[i][j].flowkey);
        if (value >= threshold)
          result[heavy[i][j].flowkey] = value;
      }
  return result;
}

template <int32_t key_len, typename hash_t>
int32_t ElasticSketch<key_len, hash_t>::compress(int32_t new_num_l_packet, bool method){
  if (n_l_packet % new_num_l_packet != 0 || new_num_l_packet != 0)
    return -1; // error
  if (n_l_packet == new_num_l_packet)
    return 0;
  int32_t (*func)(int32_t, int32_t) = (method ? sketchMax : sketchSum);
  lightpacket_t** new_mem;
  allocMem(n_l_hash, new_num_l_packet);
  for (int i = 0; i < n_l_hash; i++){
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

template <int32_t key_len, typename hash_t>
size_t ElasticSketch<key_len, hash_t>::size() const {
  return sizeof(*this)                // Instance
         + n_h_packet * sizeof(heavypacket_t<key_len>*)  // heavy
         + n_l_hash * sizeof(lightpacket_t*)    // light
         + n_h_packet * n_h_size * sizeof(heavypacket_t<key_len>) //heavy[0]
		 + n_l_hash * n_l_packet * sizeof(lightpacket_t) //light[0]
}

template <int32_t key_len, typename hash_t>
void ElasticSketch<key_len, hash_t>::clear(){
  std::fill(heavy[0], heavy[0] + n_h_packet * n_h_size * sizeof(heavypacket_t<key_len>), 0);
  std::fill(light[0], light[0] + n_l_hash * n_l_packet * sizeof(lightpacket_t), 0);
}

} // namespace OmniSketch::Sketch

#undef lightpacket_t
#undef SKETCH_INF