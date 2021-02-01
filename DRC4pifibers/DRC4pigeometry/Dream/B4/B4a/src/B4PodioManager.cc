#include "B4PodioManager.hh"
#ifdef G4MULTITHREADED
#include "G4Threading.hh"
#endif

#include <string>

B4PodioManager * B4PodioManager::m_inst_ = NULL;  

B4PodioManager::B4PodioManager():
  m_filename_prefix("simhits."),
  m_filename_suffix("podio.root")
{}

B4PodioManager * B4PodioManager::Instance()
{
  if (m_inst_ == NULL){
    m_inst_ = new B4PodioManager();
  }
  return m_inst_;
}

podio::EventStore * B4PodioManager::GetEvtStore()
{
  int threadId =  -10; // some default value different from -1, which is reserved for master in multithread mode
  
#ifdef G4MULTITHREADED
  threadId = G4Threading::G4GetThreadId();
#endif

  if (threadId == -1) {// this is the master, nothing to be done
    return NULL;
  }
  
  // Look if this thread has already a store in the map. If not, create it
  if (m_map_store.find(threadId)  == m_map_store.end()){ 

    // If it is not there, we need to create it
    podio::EventStore * l_evtstore = new podio::EventStore();
    // build the file name
    G4String filename_id = "";
    if (threadId != -10){
      filename_id = "_";
      filename_id += std::to_string(threadId);
      filename_id += "_";
    }
    G4String filename = m_filename_prefix + filename_id + m_filename_suffix;
    
    G4cout << "Podio output file name " << filename << G4endl;
    
    podio::ROOTWriter * l_writer = new podio::ROOTWriter(filename,l_evtstore);
    
    m_map_store[threadId] = l_evtstore;
    m_map_writer[threadId] = l_writer;
  }
  
  return m_map_store[threadId];
}

podio::ROOTWriter * B4PodioManager::GetWriter()
{
  G4cout << "In B4PodioManager::GetWriter()" << G4endl;
  podio::EventStore * l_eventstore = this->GetEvtStore();
  if (l_eventstore == NULL) return NULL;
  int threadId =  -10; // some default value different from -1, which is reserved for master in multithread mode
  
#ifdef G4MULTITHREADED
  threadId = G4Threading::G4GetThreadId();
#endif
  return m_map_writer[threadId];
}
  

bool B4PodioManager::Finish()
{
  G4cout << "In B4PodioManager::Finish()" << G4endl;
  
  podio::EventStore * l_eventstore = this->GetEvtStore();
  if (l_eventstore != NULL) {
    
    int threadId =  -10; // some default value different from -1, which is reserved for master in multithread mode
    
#ifdef G4MULTITHREADED
    threadId = G4Threading::G4GetThreadId();
#endif
    m_map_writer[threadId]->finish();
  }
  return true;
}
