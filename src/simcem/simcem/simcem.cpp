#include "simcem.hpp"

namespace simcem {
  constexpr const char* const Objective_t::strings[];
  constexpr const char* const Reference_t::strings[];
  constexpr const char* const Property_t::strings[];
  constexpr const char* const T_unit_t::strings[];
  constexpr const char* const Q_unit_t::strings[];
  constexpr const char* const E_unit_t::strings[];
  constexpr const char* const P_unit_t::strings[];
  constexpr const char* const L_unit_t::strings[];


  namespace detail {    
    std::tuple<std::string, std::string, std::string> tokenizeComponentID(std::string ID, std::string reqSrc, std::string reqPhase) {
      //Try parsing the ID to detect a src at the very end
      std::string sourceID = "";
      size_t src_split_idx = ID.find_first_of("@");
      if (src_split_idx != std::string::npos) {
	sourceID = ID.substr(src_split_idx+1, ID.size() - src_split_idx - 1);
	ID = ID.substr(0, src_split_idx);
      }
      
      //Try parsing the ID to detect a phase (before a src)
      std::string phaseID = "";
      size_t phase_split_idx = ID.find_first_of(":");
      if (phase_split_idx != std::string::npos) {
	phaseID = ID.substr(phase_split_idx+1, ID.size() - phase_split_idx - 1);
	ID = ID.substr(0, phase_split_idx);
      }

      if (!reqSrc.empty()) {
	if (sourceID.empty())
	  sourceID = reqSrc;
	else
	  if (sourceID != reqSrc)
	    stator_throw() << "Source name of " << reqSrc << " required, but ID requested contains alternative specification " << ID;
      }
      
      if (!reqPhase.empty()) {
	if (phaseID.empty())
	  phaseID = reqPhase;
	else
	  if (phaseID != reqPhase)
	    stator_throw() << "Phase name of " << reqPhase << " required, but ID requested contains alternative specification " << ID;
      }

      return std::make_tuple(ID, phaseID, sourceID);
    }
  }
  
  void Phase::registerIsobar(shared_ptr<Isobar> ib) {
    ib->getCurve()->setMass(_mol.mass());
    push_back(ib);
  }

  bool Bond::operator==(const Bond& b) const {
    return (((*_atom1 == *b._atom1) && (*_atom2 == *b._atom2))
	    || ((*_atom1 == *b._atom2) && (*_atom2 == *b._atom1)))
      && (_order == b._order)
      ;
  }

  void Bond::xml(Node xml) const {
    xml.add_node("Bond")
      .add_attribute("ID1", _atom1->_ID)
      .add_attribute("ID2", _atom2->_ID)
      .add_attribute("Order", _order);
  }  
}
