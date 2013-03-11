/*
 * Copyright (C) 2008-2013 CAPS entreprise.  All Rights Reserved.
 * 
 * The source code contained or described herein and all documents  related
 * to the source code ("Material") are owned  by  CAPS  entreprise  or  its
 * suppliers or licensors.
 * 
 * Title to the Material remains with CAPS entreprise or its suppliers  and
 * licensors.  The Material contains trade secrets and proprietary and con-
 * fidential information of CAPS entreprise or its suppliers and licensors.
 * 
 * The Material is protected by the French intellectual property code,  in-
 * tellectual property laws and international treaties.  No part of the Ma-
 * terial may be used, copied, reproduced, modified,  published,  uploaded,
 * posted, transmitted, distributed or disclosed in any  way  without  CAPS
 * entreprise's prior express written permission.
 * 
 * No license under any patent, copyright, trade secret or other  intellec-
 * tual property right is granted to or conferred upon you by disclosure or
 * delivery of the Material, either expressly, by implication,  inducement,
 * estoppel or otherwise.
 * 
 * Any license under such intellectual property rights  must  be  expressed
 * and approved by CAPS entreprise in writing.
 */
#ifndef OPENACCI_CONTEXT_H
#define OPENACCI_CONTEXT_H

#include <openacci/openacci.h>
#include <openacci/Region.h>
#include <openacci/Codelet.h>

#include <cstdarg>
#include <set>

namespace openacci
{
  /// Thread-local storage context
  class Context
  {
  public:
    OPENACCI_API
    Context();
  
    OPENACCI_API
    ~Context();
  
    OPENACCI_API
    static Context * getInstance();
  
    OPENACCI_API
    void init(acc_device_t);

    OPENACCI_API
    void setDeviceType(acc_device_t, int num = 0);
  
    OPENACCI_API
    acc_device_t getDeviceType() const;
  
    OPENACCI_API
    bool onDeviceType(acc_device_t dt) const;
  
    OPENACCI_API
    int getDeviceNum(acc_device_t dt) const;
  
    OPENACCI_API
    void shutdown(acc_device_t type);
  
    OPENACCI_API
    hmpprt::Device * getDevice();
  
    OPENACCI_API
    void call(const char * file_name,
              int          line_number,
              const char * function_name)
    {
      call(file_name, line_number, function_name, function_name);
    }

    // Call to an accelerated region (kernels or parallel)
    OPENACCI_API
    void call(const char * file_name,
              int          line_number,
              const char * module_name,
              const char * function_name);
  
    // Inform the runtime that we issue a fallback instead of a call (for logging only)
    OPENACCI_API
    void fallback(const char * file_name,
                  int          line_number,
                  const char * function_name);
  
    // Enter a region (declare, data, hostdata, parallel or kernels)
    OPENACCI_API
    void enterRegion(const char * file_name,
                     int          line_number,
                     int          region_kind,  // capsacc::Kind
                     int          num_arguments,
                     int          async_mode, // AsyncMode
                     int          queue_id);

    // Enter a global region (declare in a module)
    OPENACCI_API
    void enterGlobalRegion(const char * file_name,
                           int          line_number,
                           int          region_kind,  // capsacc::Kind
                           int          num_arguments,
                           int          async_mode, // AsyncMode
                           int          queue_id);
  
    // Leave a previously entered global region
    OPENACCI_API
    void leaveGlobalRegion(const char * file_name,
                           int          line_number);

    // Leave a previously entered region
    OPENACCI_API
    void leaveRegion(const char * file_name,
                     int          line_number);
  
    // Check if a module is initialized
    OPENACCI_API
    bool isInitialized(const void * module_address);

    // Record a module as initialized
    OPENACCI_API
    void recordModule(const void * module_address);

    // Add a data into the global region.
    OPENACCI_API
    void pushGlobalData(const char * file_name,
                        int          line_number,
                        const char * variable_name,
                        const void * host_address,
                        size_t       start,
                        size_t       length,
                        size_t       element_size,
                        int          transfer_mode);

    // Add a data into the current region.
    // For kernels or parallel regions, the calls must be in right order and number
    // for calling the extracted function.
    OPENACCI_API
    void pushData(const char * file_name,
                  int          line_number,
                  const char * variable_name,
                  const void * host_address,
                  size_t       start,
                  size_t       length,
                  size_t       element_size,
                  int          transfer_mode);

    // 'update' directive
    OPENACCI_API
    void updateDatas(const char   *  file_name,
                     int             line_number,
                     int             async_mode,
                     int             queue_id,
                     int             nb_variables,
                     const char   ** variable_names,
                     const void   ** host_addresses,
                     const size_t *  starts,
                     const size_t *  lengths,
                     const size_t * elements_sizes,
                     const int    *  update_sides);
 
    /// To replace a call to fortran allocate intrinsic
    OPENACCI_API
    void allocateDeviceResident(const char * file_name,
                                const int  & line_number,
                                const void * desc_address,
                                const int  & num_argument,
                                va_list      ap);
  
    OPENACCI_API
    void deallocateDeviceResident(const char * file_name,
                                  const int  & line_number,
                                  const void * desc_address);
  
    /// 'wait' directive
    OPENACCI_API
    void wait(const char * file_name,
              int          line_number,
              int          async_mode,
              int          queue_id);
  
    OPENACCI_API
    bool test(const char * file_name,
              int          line_number,
              int          async_mode,
              int          queue_id);
  
    OPENACCI_API
    void * malloc(size_t size);

    OPENACCI_API
    void   free(void *ptr);
  
    OPENACCI_API
    void * getDevicePointer(const char * file_name,
                            int          line_number,
                            const void * host_address);
    
  private:
    static void createInstance(Context ** instance);
  
    hmpprt::Device * device_;
  
    /// Available codelets
    CodeletMap codelet_map_;
    Codelet * lookupCodelet(const char * module_name, const char * function_name);
  
    // Available mirrors
    MirrorMap mirror_map_;
    MirrorMap::iterator lookupMirror(void * host_address);
  
    /// Stack of regions
    std::list<Region> region_stack_;

    /// Stack of global regions
    std::list<Region> global_region_stack_;
  
    /// Used to record buffers manually allocated through acc API (for further cleanup)
    typedef std::map<void *, hmpprt::Data *> DevicePtrMap;
    DevicePtrMap device_ptr_map_;
  
    /// User-queues
    typedef std::map<int, hmpprt::Queue *> QueueMap;
    QueueMap queue_map_;
  
    /// Anonymous queues
    typedef std::vector<hmpprt::Queue *> QueueList;
    QueueList queue_list_;

    std::set<const void*> initialized_modules_;
  
    /// \return the queue related to the given mode and id, or 0 if no queue needed.
    hmpprt::Queue * getQueue(int async_mode, int queue_id, std::string & qdesc);
  
    acc_device_t device_type_;
    int          device_num_;
  
    void clean();

    void updateData(const char * file_name,
                    int          line_number,
                    const char * variable_name,
                    const void * host_address,
                    size_t       start,
                    size_t       length,
                    size_t       element_size,
                    int          update_side,
                    int          async_mode,
                    int          queue_id);

    void enterBasicRegion(const char *        file_name,
                          int                 line_number,
                          int                 region_kind,
                          int                 num_arguments,
                          int                 async_mode,
                          int                 queue_id,
                          std::list<Region> & region_stack);

    void leaveBasicRegion(const char       *  file_name,
                          int                 line_number,
                          std::list<Region> & region_stack);


    void pushDataInRegion(const char * file_name,
                          int          line_number,
                          const char * variable_name,
                          const void * host_address_,
                          size_t       start,
                          size_t       length,
                          size_t       element_size,
                          int          transfer_mode,
                          Region     * region);

    
    /// phmpp purpose
    void makePhmppEnterLeaveRegionEvent(const char * file_name,
                                        int          line_number,
                                        int          region_kind,
                                        int          async_mode,
                                        unsigned int instance_id,
                                        bool         enter);

    void makePhmppWaitEvent(const char * file_name,
                            int          line_number,
                            int          queue_id,
                            unsigned int instance_id,
                            bool         all,
                            bool         start);

    void makePhmppEventUpdate(const char   *  file_name,
                              int             line_number,
                              int             event_start,
                              int             async_mode,
                              int             queue_id,
                              unsigned int    instance_id,
                              int             nb_variables,
                              const char   ** var_names,
                              const size_t *  starts,
                              const size_t *  lengths,
                              const size_t *  elem_sizes,
                              const int    *  update_sides);

    unsigned int newPhmppInstanceID() {static unsigned int iid = 0; iid++; return iid;}
    
  };

  
  /// class for a local Declare
  class LocalDeclare
  {
  public:
    OPENACCI_API
    LocalDeclare(const char * file_name,
                 int          line_number,
                 int          region_kind,  // capsacc::Kind
                 int          num_arguments,
                 int          async_mode, // AsyncMode
                 int          queue_id)
                 : file_name(file_name), line_number(line_number)
    {
      ctxt = Context::getInstance();
      ctxt->enterRegion(file_name,
                        line_number,
                        region_kind,
                        num_arguments,
                        async_mode,
                        queue_id);
    }
    
    OPENACCI_API
    ~LocalDeclare()
    {
      ctxt->leaveRegion(file_name,
                        line_number);
    }
    
    // Add a data into the global region.
    OPENACCI_API
    void pushGobalData(const char * file_name,
                       int          line_number,
                       const char * variable_name,
                       const void * host_address,
                       size_t       start,
                       size_t       length,
                       size_t       element_size,
                       int          transfer_mode)
    {
      ctxt->pushGlobalData(file_name,
                           line_number,
                           variable_name,
                           host_address,
                           start,
                           length,
                           element_size,
                           transfer_mode);
    }

    // Add a data into the current region.
    // For kernels or parallel regions, the calls must be in right order and number
    // for calling the extracted function.
    OPENACCI_API
    void pushData(const char * file_name,
                  int          line_number,
                  const char * variable_name,
                  const void * host_address,
                  size_t       start,
                  size_t       length,
                  size_t       element_size,
                  int          transfer_mode)
    {
      ctxt->pushData(file_name,
                     line_number,
                     variable_name,
                     host_address,
                     start,
                     length,
                     element_size,
                     transfer_mode);
    }
    
  private:
    const char * file_name;
    int          line_number;
    Context    * ctxt;
  };
} // openacci namespace


#endif // OPENACCI_CONTEXT_H

