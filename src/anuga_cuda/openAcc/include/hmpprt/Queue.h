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
#ifndef HMPPRT_QUEUE_H
#define HMPPRT_QUEUE_H

#include <hmpprt/Common.h>
#include <hmpprt/Device.h>
#include <hmpprt/Data.h>
#include <hmpprt/DeviceManager.h>
#include <hmpprt/ArgumentList.h>

namespace hmpprt
{
  class Operation
  {
  public:
    virtual ~Operation();
    virtual const char * name() const = 0;
    virtual std::string dump() const = 0;
    virtual void run(Queue * parent) = 0;
  };

  class OperationList
  {
  public:
    OperationList()
      : size_(0),
        top_(0),
        heap_()
    { }
    void * allocate();
    void push() { ++size_; }
    void clear() { size_ = 0; top_ = 0; }
    bool empty() { return size_ == 0; }
    bool pop(Operation *);

  private:
    int size_;
    int top_;
    std::vector<char> heap_;
  };

  class Queue : private NonCopyable
  {
  public:
    HMPPRT_API
    Queue();

    HMPPRT_API
    ~Queue();

    /// Launch the execution of queue operations and returns immediately.
    /// Once the queue is running, this is still possible to enqueue operations.
    HMPPRT_API
    void start();

    /// Wait until every operation in the queue is finished.
    /// Each time an operation is executed, it is removed from the queue.
    /// \param timeout_ms is in miliseconds
    /// \return false if timeout expires
    HMPPRT_API
    bool wait(size_t timeout_ms);

    /// Wait until every operation in the queue is finished.
    /// Each time an operation is executed, it is removed from the queue.
    HMPPRT_API
    void wait();

    /// Synchronous execution of the queue (equivalent to start() followed by wait(), but more efficient)
    /// Each time an operation is executed, it is removed from the queue.
    HMPPRT_API
    void execute();

    /// Clear the content of this queue, removing every operations without executing them.
    HMPPRT_API
    void clear();

    /// Return true is there is no more operation to run in the queue.
    HMPPRT_API
    bool finished();

    //@{
    /// Asynchronous operation enqueueing
    HMPPRT_API
    void enqueueAcquire(Device * d);
    HMPPRT_API
    void enqueueRelease(Device * d);
    HMPPRT_API
    void enqueueCall(Device * d, Codelet * c, ArgumentList & a);
    HMPPRT_API
    void enqueueAllocate(Data * d);
    HMPPRT_API
    void enqueueReshape(Data * dst, Data * src);
    HMPPRT_API
    void enqueueFree(Data * d);
    HMPPRT_API
    void enqueueUpload(Data * d, const void * ha);
    HMPPRT_API
    void enqueueUpload(Data * d, const void * ha, size_t offset, size_t size);
    HMPPRT_API
    void enqueueDownload(Data * d, void * ha);
    HMPPRT_API
    void enqueueDownload(Data * d, void * ha, size_t offset, size_t size);
    HMPPRT_API
    void enqueueCopy(Data * dst, Data * src);
    HMPPRT_API
    void enqueueExecute(Queue * q);
    HMPPRT_API
    void enqueueStart(Queue * q);
    HMPPRT_API
    void enqueueWait(Queue * q);
    HMPPRT_API
    void enqueueWait(Queue * q, size_t timeout_ms);
    HMPPRT_API
    void enqueueClear(Queue * q);
    HMPPRT_API
    void enqueueDelete(Data * d);
    //@}

    /// Create the thread pool (if not called, this is done at the first call to the start() method)
    static void create_thread_pool();

    /// \internal
    //@{
    void * getHandle();
    void registerHandle(Device *, void *);
    //@}

    /// \internal
    /// Run operations and awake any waiting threads after.
    /// This is an entry point for the thread running this queue.
    void runOperations();

  private:
    bool runOneOperation();
    void propagateException(void *);

    void * sem_;
    void * mutex_;
    int waiters_;
    void * handle_;
    void * result_;
    bool started_;
    Device * device_;
    OperationList operations_;
  };

  //@{
  /// Synchronous functions (for convenance only)
  inline void acquire(Device * d) { d->acquire(); }
  inline void release(Device * d) { d->release(); }
  inline void call(Device * d, Codelet * c, ArgumentList & a) { d->call(c, a); }
  inline void allocate(Data * d) { d->allocate(); }
  inline void reshape(Data * dst, Data * src) { dst->setAddress(src->getAddress()); }
  inline void free(Data * d) { d->free(); }
  inline void upload(Data * d, const void * ha) { d->upload(ha); }
  inline void upload(Data * d, const void * ha, size_t offset, size_t size) { d->upload(ha, offset, size); }
  inline void download(Data * d, void * ha) { d->download(ha); }
  inline void download(Data * d, void * ha, size_t offset, size_t size) { d->download(ha, offset, size); }
         HMPPRT_API
         void copy(Data * dst, Data * src, Queue * q = 0);
  inline void execute(Queue * q) { q->execute(); }
  inline void start(Queue * q) { q->start(); }
  inline void wait(Queue * q) { q->wait(); }
  inline void wait(Queue * q, size_t timeout_ms) { q->wait(timeout_ms); }
  inline void clear(Queue * q) { q->clear(); }
  //@}
}

#endif // HMPPRT_QUEUE_H
