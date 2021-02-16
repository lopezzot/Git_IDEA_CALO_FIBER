#ifndef B4PODIOMANAGER_H
#define B4PODIOMANAGER_H

#include <map>

#include "podio/EventStore.h"
#include "podio/ROOTWriter.h"

#include "globals.hh"

#include "G4AutoLock.hh"

namespace { G4Mutex B4PodioActionMutex = G4MUTEX_INITIALIZER; }

class B4PodioManager
{
 public:
  static B4PodioManager * Instance();
  void SetFilePrefix (G4String filename){m_filename_prefix = filename;}
  void SetFileSuffix(G4String filename){m_filename_suffix = filename;}
  bool Finish();
  podio::EventStore * GetEvtStore();
  podio::ROOTWriter * GetWriter();
    

 protected:
  G4String m_filename_prefix;
  G4String m_filename_suffix; 
  std::map<int,podio::EventStore *> m_map_store;
  std::map<int,podio::ROOTWriter *> m_map_writer;  
 private:
  static B4PodioManager * m_inst_;
  B4PodioManager();

};

#endif 
