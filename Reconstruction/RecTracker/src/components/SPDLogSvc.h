#ifndef RECTRACKER_SPDLOGSVC_H
#define RECTRACKER_SPDLOGSVC_H


// Gaudi
#include "GaudiKernel/Service.h"
#include "GaudiKernel/IService.h"
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/IMessageSvc.h"
//#include "GaudiKernel/CommonMessaging.h"
#include "GaudiKernel/MsgStream.h"

#include <iostream>
#include <memory>


#include "spdlog/spdlog.h"
#include "spdlog/details/null_mutex.h"
#include "spdlog/sinks/base_sink.h"

#include <cstdio>
#include <memory>
#include <mutex>

namespace spdlog { namespace sinks {


template <class Mutex> class gaudi_sink SPDLOG_FINAL : public base_sink<Mutex>
{
    using MyType = gaudi_sink<Mutex>;

public:
    gaudi_sink(MsgStream& aMsgStream): m_spdMess(aMsgStream)  {
    }

    static std::shared_ptr<MyType> instance()
    {
        static std::shared_ptr<MyType> instance = std::make_shared<MyType>();
        return instance;
    }

protected:
    void _sink_it(const details::log_msg &msg) override
    {
        m_spdMess << MSG::INFO << msg.raw.str() << endmsg;
    }

    void _flush() override
    {
    }
    private:
    MsgStream m_spdMess;
};
using gaudi_sink_mt = gaudi_sink<std::mutex>;
using gaudi_sink_st = gaudi_sink<details::null_mutex>;
}}

class SPDLogSvc : public extends1<Service, IService> {
public:
  /// Standard constructor
  explicit SPDLogSvc(const std::string& aName, ISvcLocator* aSL);
  /// Standard destructor
  virtual ~SPDLogSvc();
  virtual StatusCode initialize() final;
  virtual StatusCode finalize() final;

private:
};

#endif /* RECTRACKER_SPDLOGSVC_H */
