#pragma once
#include <simcem/curve.hpp>
#include <Eigen/Eigen>
#include <unordered_map>
#include <iomanip>
#include <map>
#include <sstream>
#include <fstream>
#include <vector>
#include <iostream>
#include <exception>
#include <iterator>
#include <algorithm>
#ifdef SIMCEM_NLOPT
# include <nlopt.hpp>
#endif
#define HAVE_CSTDDEF
#include <IpTNLP.hpp>
#include <IpIpoptApplication.hpp>
#include <stator/string.hpp>
#include <stator/symbolic/ad.hpp>

//Fixes for MSVC non-standardness
#include <ciso646>

/*!
  \mainpage 
  
  Simcem is written in C++ and this is documented here. The
  python interface is a simple wrapper over the C++ code, so it shares
  identical functions and class names.
*/

/*! \brief The Simcem namespace.*/
namespace simcem {
  namespace detail {    
    std::tuple<std::string, std::string, std::string> tokenizeComponentID(std::string ID, std::string reqSrc="", std::string reqPhase="");
  }
  
  class Database;
  
  /*!\brief A vector of amounts (moles or molecules) indexed by the
    species or element name (represents \f$\left\{N_i\right\}\f$).
      
    This class provides the basic molar balancing features of
    Simcem. Components can be added/subtracted from one
    another, and scaled.

    This uses a map for underlying storage as component order must
    be preserved (E.g., the equilibrium finder assigns an index to
    each variable which must not change during optimisation).
  */
  class Components : public std::map<std::string, double> {
    typedef std::map<std::string, double> Base;

  public:
    using Base::Base;

    /*! 
      @name Component arithmetic
      @{*/
    /*! \brief Arithmetic operation for manipulating the component amounts.x*/
    Components& operator+=(const Components& o) {
      for (const auto& entry: o)
	(*this)[entry.first] += entry.second;
      return *this;
    }

    Components operator+(const Components& o) const {
      return Components(*this) += o;
    }

    Components& operator-=(const Components& o) {
      for (const auto& entry: o)
	(*this)[entry.first] -= entry.second;
      return *this;
    }

    Components operator-(const Components& o) const {
      return Components(*this) -= o;
    }

    Components operator*(const double a) const {
      Components retval(*this);
      for (auto& entry: retval)
	entry.second *= a;
      return retval;
    }
  
    Components operator/(const double a) const {
      Components retval(*this);
      for (auto& entry: retval)
	entry.second /= a;
      return retval;
    }
    /*! @}*/

    
    /*! 
      @name Component access
      @{*/
    /*! \brief Access the quantity of a particular component. */
    double operator[](std::string key) const {
      auto it = find(key);
      if (it == Base::end())
	return 0;
      return it->second;
    }

    double& operator[](std::string key)
    {
      auto it = find(key);
      if (it == Base::end())
	return Base::operator[](key);
      return it->second;
    }
    /*! @} */

    
    Base::const_iterator find(std::string key) const {
      std::string ID, phase_ID, src_ID;
      std::tie(ID, phase_ID, src_ID) = detail::tokenizeComponentID(key);
      
      for (auto it=this->begin(); it != this->end(); ++it) {
	std::string trial_ID, trial_phase_ID, trial_src_ID;
	std::tie(trial_ID, trial_phase_ID, trial_src_ID) = detail::tokenizeComponentID(it->first);
	if (ID not_eq trial_ID) continue;
	if (not phase_ID.empty() and phase_ID not_eq trial_phase_ID) continue;
	if (not src_ID.empty() and src_ID not_eq trial_src_ID) continue;
	return it;
      }

      return Base::end();
    }

    Base::iterator find(std::string key) {
      std::string ID, phase_ID, src_ID;
      std::tie(ID, phase_ID, src_ID) = detail::tokenizeComponentID(key);
      
      for (auto it=this->begin(); it != this->end(); ++it) {
	std::string trial_ID, trial_phase_ID, trial_src_ID;
	std::tie(trial_ID, trial_phase_ID, trial_src_ID) = detail::tokenizeComponentID(it->first);
	if (ID not_eq trial_ID) continue;
	if (not phase_ID.empty() and phase_ID not_eq trial_phase_ID) continue;
	if (not src_ID.empty() and src_ID not_eq trial_src_ID) continue;
	return it;
      }

      return this->end();
    }
    
    /*! \brief Remove all components below a certain quantity. */
    void removeSmallComponents(const double max = 0.0) {
      for (auto it = Base::begin(); it != Base::end();)
	if (it->second <= max)
	  it = Base::erase(it);
	else
	  ++it;
    }

    /*! \brief The total amount of all components. */
    double N() const {
      double sum(0);
      for (const auto& entry: *this)
	sum += entry.second;
      return sum;
    }

    /*! \brief Returns the fractions of each component (total quantity normalised to 1). */
    Components normalised() const {
      double tot = N();
      //Prevent divide by zero
      if (tot == 0) tot += 1;
      return *this / tot;
    }

    /*! \brief Generate a string representation of the Components. */
    std::string str() const {
      std::ostringstream os;
      os << "Components{" << std::setprecision(std::numeric_limits<double>::digits10 + 1);
      for (const auto& comp : *this)
	os << "\"" << comp.first << "\":" << comp.second << ", ";
      std::string output = os.str();
      return output.substr(0, output.size()-2)+"}";
    }

    /*! \brief Return the average molar mass of the Components. */
    double m(const Database& db) const {
      return M(db) / N();
    }

    /*! \brief Determine the total mass of the Components. */
    double M(const Database& db) const;

    /*! \brief Return the elemental composition of the Components. */
    Components elements(const Database& db) const;

    /*! \brief Return the reference elemental composition of the Components. */
    Components ref_elements(const Database& db) const;
    
    static Components fromMasses(const Database& db, const Components& masses);
  protected:  
  };

  /*! \brief Output operator for Components. 
    \relates Components
  */
  inline std::ostream& operator<<(std::ostream& os, const Components& c) {
    return os << c.str();
  }
  
  /*! \brief Arithmetic operator. 
    \relates Components
  */
  inline Components operator*(const double a, const Components& o)
  { return o * a; }

  /*! \brief Data structure for an isotope of a particular element. */
  struct Isotope {
    Isotope(size_t Z, size_t N, double mass, double mass_uncertainty, double abundance, std::string symbol, std::string name, std::string category):
      _Z(Z), _N(N), _mass(mass), _mass_uncertainty(mass_uncertainty), _abundance(abundance), _symbol(symbol), _name(name), _category(category)
    {}
    
    Isotope(Node xml, std::string symbol, std::string name, std::string category):
      _Z(xml.getAttribute("Z").as<size_t>()),
      _N(xml.getAttribute("N").as<size_t>()),
      _mass(xml.getAttribute("Mass").as<double>() * 1e-3),
      _mass_uncertainty(xml.getAttribute("MassUncertainty").as<double>() * 1e-3),
      _abundance(xml.getAttribute("P").as<double>()),
      _symbol(symbol),
      _name(name),
      _category(category)
    {}

    void xml(stator::xml::Node node) const {
      node
	.add_attribute("Z", _Z)
	.add_attribute("N", _N)
	.add_attribute("Mass", _mass * 1e3)
	.add_attribute("MassUncertainty", _mass_uncertainty * 1e3)
	.add_attribute("P", _abundance)
	;
    }

    /*! \brief Number of protons in the element. */
    size_t _Z;
    /*! \brief Neutron count. */
    size_t _N;
    /*! \brief Mass. */
    double _mass;
    /*! \brief Uncertainty in the mass value. */
    double _mass_uncertainty;
    /*! \brief Natural abundance of the isotope (where available). */
    double _abundance;
    /*! \brief Short symbol of the element (e.g., Cl for Chlorine). */
    std::string _symbol;
    /*! \brief Full name of the element. */
    std::string _name;
    /*! \brief Category of the element. */
    std::string _category;
  };
    
  /*! \brief Data structure for a single element. */
  struct Element : public std::vector<Isotope>, Isotope {
    typedef std::vector<Isotope> Base;
    
    Element(size_t Z, size_t N, double mass, double mass_uncertainty, double abundance, std::string symbol, std::string name, std::string category, Base isotopes, size_t group, size_t period, std::string block="", std::string referenceComponentID=""):
      Isotope(Z, N, mass, mass_uncertainty, abundance, symbol, name, category),
      Base(isotopes),
      _group(group),
      _period(period),
      _block(block),
      _referenceComponentID(referenceComponentID)
    {}
    
    Element(Node xml):
      Isotope(xml, xml.getAttribute("Symbol"), xml.getAttribute("Name"), xml.getAttribute("Category")),
      _group(xml.getAttribute("Group").as<size_t>()),
      _period(xml.getAttribute("Period").as<size_t>())
    {
      if (xml.hasAttribute("Block"))
	_block = xml.getAttribute("Block");
      for (Node isonode = xml.findNode("Isotope"); isonode.valid(); ++isonode)
	this->push_back(Element::Isotope(isonode, _symbol, _name, _category));

      if (xml.hasNode("ReferenceComponent"))
	_referenceComponentID = xml.getNode("ReferenceComponent").getAttribute("ID");
    }

    void xml(stator::xml::Node node) const {
      stator::xml::Node xml_element = node.add_node("Element");

      Isotope::xml(xml_element);

      xml_element
	.add_attribute("Name", _name)
	.add_attribute("Symbol", _symbol)
	.add_attribute("Category", _category)
	.add_attribute("Group", _group)
	.add_attribute("Period", _period)
	.add_attribute("Block", _block)
	;

      if (!_referenceComponentID.empty())
	xml_element.add_node("ReferenceComponent").add_attribute("ID", _referenceComponentID);
      
      for (const auto& isotope: *this){
	stator::xml::Node xml_isotope = xml_element.add_node("Isotope");
	isotope.xml(xml_isotope);
      }
    }

    /*! 
      @name Periodic table data
      @{*/
    /*! \brief Periodic table group. */
    size_t _group;
    /*! \brief Periodic table period. */
    size_t _period;
    /*! \brief Periodic table block. */
    std::string _block;
    /*! @} */
    /*! \brief Periodic table block. */
    std::string _referenceComponentID;
  };      
  
  struct Isobar {
    Isobar(shared_ptr<Curve> curve, double p):
      security("Public"),
      _curve(curve),
      _p(p)
    {}
      
    Isobar(Node xml):
      security("Public"),
      _curve(Curve::load(xml, Objective_t::T)),
      _p(xml.getAttribute("p").as<double>())
    {
      
      if (xml.hasAttribute("Security"))
	security = xml.getAttribute("Security").as<std::string>();
    }

    virtual ~Isobar() {}
    
    void xml(stator::xml::Node node) const {
      stator::xml::Node xml_isobar = node.add_node("Isobar");
      _curve->xml(xml_isobar);
      xml_isobar.add_attribute("p", _p);
      if (security != "Public")
	xml_isobar.add_attribute("Security", security);
    }
    
    double p() const {
      return _curve->normalise_units(Objective_t(Objective_t::p), _p);
    }

    shared_ptr<Curve> getCurve() {
      return _curve;
    }

    void setCurve(shared_ptr<Curve> c) {
      _curve = c;
    }

    std::string security;
    
  protected:
    shared_ptr<Curve> _curve;
    double _p;
  };

  struct Component;

  /*! \brief Data structure for a phase of a particular molecule. */
  struct Phase : public std::vector<shared_ptr<Isobar> > {
    void load_xml(Node xml) {
      for (Node inode = xml.findNode("Alias"); inode.valid(); ++inode)
	_aliases.push_back(inode.getAttribute("Name"));
      
      for (Node inode = xml.findNode("Isobar"); inode.valid(); ++inode)
	registerIsobar(make_shared<Isobar>(inode));

      if (xml.hasNode("Comments"))
	_comments = xml.getNode("Comments").getValue();
    }
    
    Phase(std::string name, std::string type, Component& md):
      _name(name), _type(type),
      _mol(md)
    {}

    void xml(stator::xml::Node node) const {
      stator::xml::Node xml_phase = node.add_node("Phase");
      for (const auto& alias: _aliases)
	xml_phase.add_node("Alias").add_attribute("Name", alias);
      
      xml_phase
	.add_attribute("Name", _name)
	.add_attribute("Type", _type)
	;

      if (!_comments.empty())
	xml_phase.add_node("Comments") = _comments;
      
      for (const auto& isobar: *this)
	isobar->xml(xml_phase);
    }

    void registerIsobar(shared_ptr<Isobar> ib);

    void add_alias(std::string alias) {
      std::string lwr_alias = alias;
      boost::algorithm::to_lower(lwr_alias);

      std::string lwr_name = _name;
      boost::algorithm::to_lower(lwr_name);
      if (lwr_name == lwr_alias) return;
      
      for (std::string a : _aliases) {
	boost::algorithm::to_lower(a);
	if (lwr_alias == a)
	  return;
      }
      _aliases.emplace_back(alias);
    }

    void remove_alias(std::string alias) {
      _aliases.erase(std::remove_if(_aliases.begin(), _aliases.end(), [alias](std::string a){ return boost::iequals(alias,a); }), _aliases.end());
    }
    
    const std::vector<std::string>& getAliases() const {
      return _aliases;
    }

    std::string getComments() const { return _comments; }
    void setComments(std::string comments) { _comments = comments; }
    
    std::string _name;
    std::string _type;
    std::string _comments;
    std::vector<std::string> _aliases;
    Component& _mol;
  };

  struct Atom;
  
  struct Bond {
    Bond(shared_ptr<Atom> atom1, shared_ptr<Atom> atom2, double order):
      _atom1(atom1), _atom2(atom2), _order(order)
    {}

    Bond() {}

    bool operator==(const Bond& b) const;
    
    void xml(Node xml) const;
    
    shared_ptr<Atom> _atom1;
    shared_ptr<Atom> _atom2;
    double _order;
  };
  
  struct Atom {
    Atom():
      _ID(0), _Z(0), _N(0), _charge(0), _x(0), _y(0)
    {}
    
    Atom(size_t ID, size_t Z, int N, double x, double y, double charge):
      _ID(ID), _Z(Z), _N(N), _charge(charge), _x(x), _y(y)
    {}

    Atom(Node xml):
      _ID(xml.getAttribute("ID").as<size_t>()),
      _Z(xml.getAttribute("Z").as<size_t>()),
      _x(xml.getAttribute("x").as<double>()),
      _y(xml.getAttribute("y").as<double>())
    {
      _N = xml.hasAttribute("N") ? int(xml.getAttribute("N").as<unsigned int>()) : -1;
      _charge = xml.hasAttribute("Charge") ? xml.getAttribute("Charge").as<double>() : 0;
    }
    
    void xml(Node xml) const {
      Node atomxml = xml
	.add_node("Atom")
	.add_attribute("ID", _ID)
	.add_attribute("Z", _Z)
	.add_attribute("x", _x)
	.add_attribute("y", _y)
	;
      
      if (_N != -1)
	atomxml.add_attribute("N", _N);
      
      if (_charge != 0)
	atomxml.add_attribute("Charge", _charge);
    }

    bool operator==(const Atom& a) const {
      return (_ID == a._ID)
	&& (_Z == a._Z)
	&& (_N == a._N)
	&& (_charge == a._charge)
	&& (_bonds == a._bonds)
	;
    }
    
    size_t _ID;
    size_t _Z;
    int _N;
    double _charge;
    std::vector<Bond> _bonds;
    double _x, _y;
  };
  
  /*! \brief Data structure for a single molecule. */
  struct Component : public std::vector<Phase> {
    void xml(stator::xml::Node node) const {
      stator::xml::Node xml_molecule = node.add_node("Component");

      //Required
      xml_molecule.add_attribute("Identifier", _identifier);

      if (!_formula.empty())
	xml_molecule.add_attribute("Formula", _formula);
      
      for (const auto& alias: _aliases)
	xml_molecule.add_node("Alias").add_attribute("Name", alias);
	
      for (const auto& element: _elements)
	xml_molecule
	  .add_node("Atom")
	  .add_attribute("Name", element.first)
	  .add_attribute("Quantity", element.second)
	  ;

      if (!_structure.empty()) {
	auto struct_node = xml_molecule.add_node("Structure");
	for (const auto& atom : _structure)
	  atom->xml(struct_node);

	for (const auto& atom : _structure)
	  for (const auto& bond : atom->_bonds)
	    if (bond._atom1->_ID < bond._atom2->_ID)
	      bond.xml(struct_node);
      }
      
      for (const auto& property: _properties)
	for (const auto & p : property.second)
	  p.xml(xml_molecule.add_node(property.first));
	
      for (const auto& phase: *this)
	phase.xml(xml_molecule);

    }

    Phase& registerPhase(std::string phaseID, std::string type) {
      for (auto& phase: *this) {
	if (phaseID == phase._name) {
	  if (type == phase._type) {
	    return phase;
	  } else
	    stator_throw() << "Registering phase " << phaseID << " of type " << type <<  " for Component " << _identifier << "; however, this phase already exists with type=" << phase._type;
	}
      }

      push_back(Phase(phaseID, type, *this));
      return back();
    }

    const Phase& getPhase(std::string phaseID, std::string type) const {
      for (const auto& phase: *this)
	if (type.empty() || (type == phase._type))
	  if (phaseID.empty() || (phaseID == phase._name))
	    return phase;
      stator_throw() << "Failed to find phase type=\"" << type << "\" with ID=\"" << phaseID << "\" in molecule ID=\"" << _identifier << "\"";
    }
    
    double mass() const { return _mass; }
    std::string getID() const { return _identifier; }
    const Components& getElements() const { return _elements; }
    const std::string& getFormula() const { return _formula; }
    void setFormula(std::string fm) { _formula = fm; }

    struct Property : Data {
      Property(Node xml):
	Data(xml),
	_value(xml.getAttribute("Value").as<double>()),
	_type(xml.getName())
      {}
	
      Property(Property_t type, double v, std::string source, T_unit_t Tunit = T_unit_t::K, Q_unit_t Qunit = Q_unit_t::mol, E_unit_t Eunit = E_unit_t::J, P_unit_t Punit = P_unit_t::Pa, L_unit_t Lunit = L_unit_t::m):
	Data("", source, Tunit, Qunit, Eunit, Punit, Lunit),
	_value(v),
	_type(type)
      {}

      void xml(Node xml) const {
	xml.add_attribute("Value", _value);
	Data::xml(xml);
      }

      double getValue() const {
	return normalise_units(_type, _value);
      }

      double getOrigValue() const {
	return _value;
      }
      
      Property_t getType() const {
	return _type;
      }

      bool operator==(const Property& p) const {
	return (_value == p._value)
	  && (_type == p._type)
	  && Data::operator==(p);
      }

    protected:
      double _value;
      Property_t _type;
    };

    typedef std::map<Property_t, std::vector<Property> > PropertyMap;
    
    PropertyMap& getProperties() { return _properties; }
    const PropertyMap& getProperties() const { return _properties; }

    void registerProperty(Property_t type, std::string source, double val, T_unit_t Tunit = T_unit_t::K, Q_unit_t Qunit = Q_unit_t::mol, E_unit_t Eunit = E_unit_t::J, P_unit_t Punit = P_unit_t::Pa, L_unit_t Lunit = L_unit_t::m)
    {
      _properties[type].push_back(Property(type, val, source, Tunit, Qunit, Eunit, Punit, Lunit));
    }

    
    Component(std::string ID, Components elements, const Database& db, std::string formula = ""):
      _elements(elements),
      _formula(formula),
      _identifier(ID)
    {
      if (_formula.empty())
	//Build a formula (assuming we have an elemental description,
	//if not, an empty formula is fine
	for (const auto& val: _elements)
	  _formula = _formula + val.first + boost::lexical_cast<std::string>(val.second);
      
      init_mass(db);
    }

    inline
    void load(Node xml) {
      for (Node inode = xml.findNode("Alias"); inode.valid(); ++inode)
	add_alias(inode.getAttribute("Name"));
      
      for (Node pnode = xml.findNode("Phase"); pnode.valid(); ++pnode)
	registerPhase(pnode.getAttribute("Name"), pnode.getAttribute("Type")).load_xml(pnode);
	  
      for (auto prop : Property_t::strings)
	for (Node anode = xml.findNode(prop); anode.valid(); ++anode)
	  getProperties()[Property_t(prop)].emplace_back(anode);

      if (xml.hasNode("Structure")) {
	std::vector<shared_ptr<Atom> > newstruct;
	Node struct_xml = xml.getNode("Structure");
	for (Node anode = struct_xml.findNode("Atom"); anode.valid(); ++anode)
	  newstruct.emplace_back(new Atom(anode));

	for (Node bnode = struct_xml.findNode("Bond"); bnode.valid(); ++bnode) {

	  shared_ptr<Atom> a1;
	  shared_ptr<Atom> a2;
	  for (auto& atom : newstruct) {
	    if (atom->_ID == bnode.getAttribute("ID1").as<size_t>())
	      a1 = atom;
	    if (atom->_ID == bnode.getAttribute("ID2").as<size_t>())
	      a2 = atom;
	  }
	  if (!a1 || !a2)
	    stator_throw() << "Could not find one of the Atom IDs requested.\n" << bnode.getPath();

	  a1->_bonds.push_back(Bond(a1, a2, bnode.getAttribute("Order").as<double>()));
	  a2->_bonds.push_back(Bond(a2, a1, bnode.getAttribute("Order").as<double>()));
	}
	
	if (!_structure.empty() && _structure != newstruct)
	  stator_throw() << "New structure for molecule does not match existing structure!\n"
			 << struct_xml.getPath();
	_structure = newstruct;
      }
    }
    
    void add_alias(std::string alias) {
      std::string lwr_alias = alias;
      boost::algorithm::to_lower(lwr_alias);

      std::string lwr_name = _identifier;
      boost::algorithm::to_lower(lwr_name);
      if (lwr_name == lwr_alias) return;
      
      for (std::string a : _aliases) {
	boost::algorithm::to_lower(a);
	if (lwr_alias == a)
	  return;
      }
      _aliases.emplace_back(alias);
    }


    void remove_alias(std::string alias) {
      _aliases.erase(std::remove_if(_aliases.begin(), _aliases.end(), [alias](std::string a){ return boost::iequals(alias,a); }), _aliases.end());
    }

    const std::vector<std::string>& getAliases() const {
      return _aliases;
    }

    
    const std::vector<shared_ptr<Atom> >& getStructure() const { return _structure; }
    std::vector<shared_ptr<Atom> >& getStructure()  { return _structure; }
    
  protected:
    void init_mass(const Database& db);
      
    /*! \brief The elements and their amounts which comprise the molecule.*/
    Components _elements;
    /*! \brief A representative chemical formula for the molecule. */
    std::string _formula;
    /*! \brief A unique identifier for the molecule. */
    std::string _identifier;
    /*! \brief Average mass of the molecule. */
    double _mass;
    /*! \brief Additional properties of the molecule. */
    PropertyMap _properties;
    
    std::vector<std::string> _aliases;
    std::vector<shared_ptr<Atom> > _structure;
  };

  /*! \brief A database of thermophysical data required for the
    physical models. */
  class Database
  {
  protected:
    Database() {}
    
    std::unordered_map<std::string, Element> _elements;
    std::unordered_map<std::string, Component> _components;
    
  public:
    static constexpr double pi = 3.1415926535897932384626433832795029L;

    //Taken from CODATA 2014 https://physics.nist.gov/cuu/Constants/
    double avogadro = 6.022140857e23;
    double kB = 1.38064852e-23;
    double steffanBoltzmannConstant = 5.670367e-8;
    //! \brief Standard gravity acceleration.
    constexpr static double g = 9.80665; 
    
    //Derived constants Calculated in the constructor
    double R() const { return kB * avogadro; }

    static inline Components getDryAir() {
      //Dry air, as reported on wikipedia for 2015
      return {{"N2",78.084}, {"O2",20.946},  {"Ar",0.934}, {"CO2",0.04}, {"Ne",0.001818}, {"He",0.000524}, {"CH4",0.000179}, {"Kr",0.000114}};
    }

    
    Component& registerComponent(std::string ID, Components elements, std::string formula = "") {
      auto it = _components.find(ID);
      if (it == _components.end())
	return _components.insert(std::make_pair(ID, Component(ID, elements, *this, formula))).first->second;
      else {
	//Check for consistency
	if (it->second.getElements() != elements)
	  stator_throw() << "Failed to register molecule " << ID << " as it already exists and has a different elemental composition!";
	if (!formula.empty() && (it->second.getFormula() != formula))
	  stator_throw() << "Failed to register molecule " << ID << " as its formula (" << formula << ") does not match existing formula ("<< it->second.getFormula() << ")!";
	return it->second;
      }
    }
    
    const Component& getComponent(std::string molID) const {
      std::tie(molID, std::ignore, std::ignore) = detail::tokenizeComponentID(molID);
      auto it = _components.find(molID);
      if (it == _components.end())
	stator_throw() << "Failed to find molecule " << molID;
      return it->second;
    }

    std::string xml() const {
      std::ostringstream os;
      stator::xml::Document doc;

      auto root = doc.add_node("SimCem");

      auto xml_constants = root.add_node("Constants");

      xml_constants.add_node("Avogadro").add_attribute("Value", avogadro);
      xml_constants.add_node("kB").add_attribute("Value", kB);
      xml_constants.add_node("steffanBoltzmann").add_attribute("Value", steffanBoltzmannConstant);
      
      auto xml_elements = root.add_node("Elements");      
      for (const auto& element: _elements)
	element.second.xml(xml_elements);

      auto xml_components = root.add_node("Components");
      for (const auto& molecule: _components)
	molecule.second.xml(xml_components);

      os << doc;
      return os.str();
    }

    static shared_ptr<Database> create() {
      return shared_ptr<Database>(new Database());
    }

    static shared_ptr<Database> create_from_file(std::string file) {
      shared_ptr<Database> db(new Database());
      db->load(file);
      return db;
    }

    inline
    void load(std::string filename) {
      stator::xml::Document doc(filename);
      Node mainNode = doc.getNode("SimCem");

      if (mainNode.hasNode("Constants")) {
	Node constants = mainNode.getNode("Constants");

	if (constants.hasNode("Avogadro"))
	  avogadro = constants.getNode("Avogadro").getAttribute("Value").as<double>();

	if (constants.hasNode("kB"))
	  kB = constants.getNode("kB").getAttribute("Value").as<double>();

	if (constants.hasNode("steffanBoltzmann"))
	  steffanBoltzmannConstant = constants.getNode("steffanBoltzmann").getAttribute("Value").as<double>();
      }
      
      if (mainNode.hasNode("Elements"))
	for (Node elenode = mainNode.getNode("Elements").findNode("Element"); elenode.valid(); ++elenode)
	  _elements.insert(std::unordered_map<std::string, Element>::value_type(elenode.getAttribute("Symbol"), Element(elenode)));

      if (mainNode.hasNode("Components"))
	for (Node molnode = mainNode.getNode("Components").findNode("Component"); molnode.valid(); ++molnode) {
	  Components elements;
	  for (Node anode = molnode.findNode("Atom"); anode.valid(); ++anode)
	    elements[anode.getAttribute("Name")] = anode.getAttribute("Quantity").as<double>();

	  std::string formula = "";
	  if (molnode.hasAttribute("Formula"))
	    formula = molnode.getAttribute("Formula");
	  auto & mol = registerComponent(molnode.getAttribute("Identifier"), elements, formula);

	  mol.load(molnode);
	}
    }

    const Element& getElement(std::string symbol) const {
      auto it = _elements.find(symbol);
      if (it == _elements.end())
	stator_throw() << "Failed to find element " << symbol;
      return it->second;
    }

    const std::unordered_map<std::string, Element>& getElements() const {
      return _elements;
    }

    const std::unordered_map<std::string, Component>& getComponents() const {
      return _components;
    }
  };

  inline double Components::M(const Database& db) const {
    double sum(0);
    for (const auto& entry: *this)
      sum += entry.second * db.getComponent(entry.first).mass();
    return sum;
  }
  
  inline Components Components::elements(const Database& db) const {
    Components retval;
    for (const auto& mol : *this)
      retval += mol.second * db.getComponent(mol.first).getElements();
    return retval;
  }

  inline Components Components::ref_elements(const Database& db) const {
    Components retval;
    for (const auto& mol : *this) {
      for (const auto& element : db.getComponent(mol.first).getElements()) {
	//Grab the component ID for the element reference
	const std::string ref_ele_ID = db.getElement(element.first)._referenceComponentID;
	//Figure out how many of that element are in the reference, compared to this
	const double scale = element.second / db.getComponent(ref_ele_ID).getElements()[element.first];
	retval[ref_ele_ID] += scale;
      }
    }
    return retval;
  }
  
  inline Components Components::fromMasses(const Database& db, const Components& masses) {
    Components retval;
    for (const auto& mol : masses)
      retval[mol.first] = mol.second / db.getComponent(mol.first).mass();
    return retval;
  }

  
  inline void Component::init_mass(const Database& db) {
    _mass = 0;
    //We can't use total_mass(), as this requires a look up of
    //Component, which may not be loaded yet.
    for (const auto& element : _elements)
      _mass += db.getElement(element.first)._mass * element.second;
  }

  struct IsobarIterator {
    typedef Component::const_iterator outer_iterator;
    typedef Phase::const_iterator inner_iterator;

    typedef typename inner_iterator::value_type      value_type;
    typedef typename inner_iterator::difference_type difference_type;
    typedef typename inner_iterator::pointer         pointer;
    typedef typename inner_iterator::reference       reference;
      
    outer_iterator _out_it;
    outer_iterator _out_end;
    inner_iterator _in_it;

    IsobarIterator() {}
    
    IsobarIterator(const outer_iterator it):
      _out_it(it),
      _out_end(it)
    {}

    void initialise(const Component& mol)
    {
      _out_it = mol.begin();
      _out_end = mol.end();

	if (_out_it != mol.end())
	_in_it = _out_it->begin();
      
      advance();
    }

    reference operator*()  const { return *_in_it;  }
    pointer   operator->() const { return &*_in_it; }
      
    IsobarIterator& operator++() {
      ++_in_it;
      advance();
      return *this;
    }

    IsobarIterator operator++(int)
    {
      IsobarIterator it(*this);
      ++*this;
      return it;
    }
      
    void advance() {
      //Skip any inner items
      while ((_in_it != _out_it->end()) && skip_inner())
	++_in_it;

      //Check if we need to skip the outer
      while ((_out_it != _out_end) && ((_in_it == _out_it->end()) || skip_outer()))
        {
	  ++_out_it;
	  if (_out_it != _out_end) {
	    _in_it = _out_it->begin();

	    //Skip required inner items
	    while ((_in_it != _out_it->end()) && skip_inner())
	      ++_in_it;
	  }
        }
    }

    friend bool operator==(const IsobarIterator& a, 
                           const IsobarIterator& b)
    {
        if (a._out_it != b._out_it)
            return false;

        if (a._out_it != a._out_end && 
            b._out_it != b._out_end &&
            a._in_it != b._in_it)
            return false;

        return true;
    }

    friend bool operator!=(const IsobarIterator& a,
                           const IsobarIterator& b)
    {
        return !(a == b);
    }

    virtual bool skip_outer() {
      return false;
    }

    virtual bool skip_inner() {
      return false;
    }
  };
  
  /*! \brief A container of Isobar instances which is filtered out of the database.
   */
  struct IsobarRange : std::vector<shared_ptr<Isobar> > {
    typedef std::vector<std::function<bool(const Phase&)>> PhaseFilters;
    typedef std::vector<std::function<bool(const shared_ptr<Isobar>&)>> IsobarFilters;
    
    IsobarRange(const Database& db, std::string ID, PhaseFilters&& pfilters={}, IsobarFilters&& ifilters={}):
      _phaseFilters(pfilters),
      _isobarFilters(ifilters)
    {
      std::string phaseID, sourceID;
      std::tie(ID, phaseID, sourceID) = detail::tokenizeComponentID(ID, sourceID, phaseID);
      const Component& mol = db.getComponent(ID);

      if (!phaseID.empty())
	_phaseFilters.push_back([=](const Phase& phase) { return phaseID != phase._name; });
      if (!sourceID.empty())
	_isobarFilters.push_back([=](const shared_ptr<Isobar>& isobar) { return isobar->getCurve()->_source != sourceID; });
	  	
      //First pass to count the number of isobars
      size_t isobar_count = 0;
      forEachIsobar(mol, [&](const shared_ptr<Isobar>&) {++isobar_count;});
      //Allocate at once, to minimise dynamic allocations
      this->resize(isobar_count);

      //Now copy out the matching isobars
      isobar_count = 0;
      forEachIsobar(mol, [&](const shared_ptr<Isobar>& isobar) {(*this)[isobar_count++] = isobar;});
	
    }

    void forEachIsobar(const Component& mol, std::function<void(const shared_ptr<Isobar>&)> func) {
      for (const Phase& phase: mol) {
	bool skip = false;
	for (const auto& filter : _phaseFilters)
	  if (filter(phase)) { skip = true; break; }
	if (skip) continue;
	  
	for (const shared_ptr<Isobar>& isobar : phase)
	  {
	    bool skip = false;
	    for (const auto& filter : _isobarFilters)
	      if (filter(isobar)) { skip = true; break; }
	    if (skip) continue;
	    func(isobar);
	  }
      }	
    }
    
    PhaseFilters _phaseFilters;
    IsobarFilters _isobarFilters;
  };
  
  /*! \brief The base class for any thermodynamic phase model.
   */
  class Model : public Components {
  public:
    Model(shared_ptr<Database> db, Components components, const std::string type):
      Components(components),
      _db(db), _type(type)
    {}

    virtual Objective_t Y(size_t idx) const = 0;

    virtual double p() const = 0;
    virtual double T() const = 0;
    
    //Extensive variables
    /*! \brief Partial molar volume */
    virtual double v(std::string molID) const = 0;

    /*! \brief Volume. */
    virtual double V() const = 0;

    /*! \brief Partial molar entropy. */
    virtual double s(std::string molID) const = 0;

    /*! \brief Entropy. */
    double S() const { return extensiveProperty<&Model::s>(); }

    std::string type() const { return _type; }
    
    /*! \brief Partial molar Gibb's free energy/chemical potential */
    virtual double chemPot(std::string molID) const = 0;

    /*! \brief Gibbs free energy. */
    double G() const { return extensiveProperty<&Model::chemPot>(); }

    /*! \brief Partial molar enthalpy. */
    virtual double h(std::string molID) const = 0;

    /*! \brief Enthalpy. */
    double H() const { return extensiveProperty<&Model::h>(); }

    /*! \brief Partial molar internal energy. */
    virtual double u(std::string molID) const = 0;
      
    /*! \brief Internal energy. */
    double U() const { return extensiveProperty<&Model::u>(); }

    /*! \brief Partial molar Helmholtz free energy. */
    virtual double a(std::string molID) const = 0;

    /*! \brief Helmholtz free energy. */
    double A() const { return extensiveProperty<&Model::a>(); }

    //Material derivatives
    /*! \brief Isobaric heat capacity. */
    virtual double Cp() const = 0;

    /*! \brief Thermal expansion coefficient. */
    virtual double Alpha() const = 0;

    /*! \brief Isothermal compressibility. */
    virtual double Beta() const = 0;

    /*! \brief Output a text string representation of this phase.
     */
    virtual std::string str() const = 0;
      
    typedef double (Model::*SpecificPropPtr_t)(std::string) const;

    double f(Objective_t f) const {
      switch (f) {
      case Objective_t::G:
	return G();
      case Objective_t::H:
	return H();
      case Objective_t::negS:
	return -S();
      case Objective_t::S:
	return S();
      case Objective_t::U:
	return U();
      case Objective_t::A:
	return A();
      case Objective_t::p:
	return p();
      case Objective_t::V:
	return V();
      default:
	stator_throw() << "Unhandled Objective_t " << f;
      }
    }

    virtual double dfdT(Objective_t) const = 0;
    virtual double dfdY2(Objective_t) const = 0;
    virtual double dfdNi(Objective_t, std::string molID) const = 0;

    shared_ptr<Database> db() const {
      return _db;
    }

    double M() const { return Components::M(*_db); }
    double m() const { return Components::m(*_db); }

    Components elements() const { return Components::elements(*_db); }

    virtual void set(Objective_t val, double target, Objective_t constant) = 0;
    
    bool hasData() const {
      for (const auto& comp : *this) {
	bool foundData = false;
	for (IsoIt it(*db(), comp.first, _type); it != it.end(); ++it)
	  if ((*it)->getCurve()->inRange(T()) && !dynamic_cast<const TabulatedCurve*>((*it)->getCurve().get())) {
	    foundData = true;
	    break;
	  }
	if (!foundData)
	  return false;
      }
      return true;
    }

  protected:      
    friend class System;
    
    struct IsoIt : IsobarIterator {
      std::string _phaseID, _sourceID, _ID, _type;
      const Database& _db;
      IsoIt(const Database& db, const std::string ID, const std::string type):
	_type(type),
	_db(db)
      {
	std::tie(_ID, _phaseID, _sourceID) = detail::tokenizeComponentID(ID);
	initialise(db.getComponent(ID));
      }

      IsobarIterator end() {
	return IsobarIterator(_db.getComponent(_ID).end());
      }
      
      virtual bool skip_outer() {
	return (_out_it->_type != _type)
	  || (!_phaseID.empty() && (_phaseID != _out_it->_name));
      }

      virtual bool skip_inner() {
	return ((*_in_it)->getCurve()->getVariable() != Objective_t::G)
	  || (!_sourceID.empty() && (_sourceID != (*_in_it)->getCurve()->_source));
      }
    };

    struct DataPoint {
      double x;
      double f;
      double df;
      double ddf;
      double p;
    };
    
    DataPoint getGdata(std::string ID, double T) const {
      shared_ptr<Isobar> iso;
      double Tdist = HUGE_VAL;
      
      for (IsoIt it(*db(), ID, _type); it != it.end(); ++it) {
	if (!dynamic_cast<const FunctionCurve*>((*it)->getCurve().get()))
	  continue;
	
	const auto& isobar = (*it);
	const auto& curve = isobar->getCurve();
	double dist = (T > curve->xmax()) * (T - curve->xmax()) + (T < curve->xmin()) * (T - curve->xmin());
	if (std::abs(dist) < std::abs(Tdist)) {
	  iso = isobar;
	  Tdist = dist;
	}
      }

      if (!iso)
	stator_throw() << "Failed to find any suitable IdealGas G isobar data @T=" << T  << " for " << ID;

	sym::Var<sym::vidx<'T'>> varT;
      if (Tdist == 0) {
	//Just evaluate the function at the requested point
	auto eval = sym::ad<2>(iso->getCurve()->getFunction(), varT = T);
	return DataPoint{T, eval[0], eval[1], eval[2] * 2, iso->p()};
      } else {
	//Evaluate by extrapolating from the bounds using a linear Cp
	double Tlim = (Tdist>0) ? iso->getCurve()->xmax() : iso->getCurve()->xmin();
	auto eval = sym::ad<3>(iso->getCurve()->getFunction(), varT = Tlim);
	
	const double cp0 = Tlim * Tlim * eval[3] * 6;
	const double cp1 = -(eval[2]*2 + cp0 / Tlim);
	const double c1 = eval[1] + (cp0 * std::log(Tlim) + cp1 * Tlim);
	const double c2 = eval[0] + cp0 * (Tlim * std::log(Tlim) - Tlim) + cp1 * Tlim * Tlim / 2 - c1 * Tlim;
	auto mu = -cp0 * (varT * sym::log(varT) - varT) - cp1 * varT * varT / 2 + c1 * varT + c2;
	auto eval2 = sym::ad<2>(mu, varT = T);
	return DataPoint{T, eval2[0], eval2[1], eval2[2] * 2, iso->p()};
      }
    }

    double cp_ig0(std::string molID) const {
      const auto data = getGdata(molID, T());
      return - T() * data.ddf;
    }
    
    template<SpecificPropPtr_t specificProp>
    double extensiveProperty() const {
      double sum(0);
      for (const auto& comp : *this)
	sum += comp.second * (this->*specificProp)(comp.first);
      return sum;
    }    
    
    shared_ptr<Database> _db;
    std::string _type;
  };


  /*! \brief A class which converts an existing model into a (infinite) .
  */
  class ModelExcess : public Model {
    shared_ptr<Model> _base_model;
    
  public:
    ModelExcess(shared_ptr<Model> model):
      Model(model->db(), Components(*model), model->type()),
      _base_model(model)
    {}
    
#define fail() stator_throw() << "Cannot use an excess model this way!"; 
    
    virtual double T() const {return _base_model->T(); }
    virtual double p() const { return _base_model->p(); }
    virtual double V() const { fail(); }
    void set(Objective_t val, double target, Objective_t constant) { fail(); }
    virtual double Cp() const { fail(); }
    virtual double Alpha() const { fail(); }
    virtual double Beta() const { fail(); }
    virtual double chemPot(std::string molID) const { fail(); }
    virtual std::string str() const { fail(); }
    virtual double s(std::string molID) const { fail(); }
    virtual double v(std::string molID) const { fail(); }
    virtual double h(std::string molID) const { fail(); }
    virtual double u(std::string molID) const { fail(); }
    virtual double a(std::string molID) const { fail(); }
    virtual Objective_t Y(size_t i) const { fail(); }
    virtual double dfdT(Objective_t f) const { fail(); }
    virtual double dfdY2(Objective_t f) const { fail(); }
    virtual double dfdNi(Objective_t p, std::string molID) const { fail(); }
  };

  /*! \brief A base class for phases using temperature as a variable.
  */
  class ModelT : public Model {
  public:
    ModelT(shared_ptr<Database> db, Components components, const double T, const std::string type):
      Model(db, components, type),
      _T(T)
    {}

    virtual double T() const { return _T; }
  protected:
    double _T;
  };
  
  inline std::ostream& operator<<(std::ostream& os, const Model& c) {
    return os << c.str();
  }

  /*! \brief A base class for phases described using
    \f$\vec{X}=\left\{T,p,\left\{N_i\right\}^{N_c}\right\}\f$.
  */
  class ModelTp : public ModelT {
  public:
    ModelTp(shared_ptr<Database> db, Components components, const double T, const double p, const std::string type):
      ModelT(db, components, T, type),
      _p(p)
    {}

    virtual double p() const { return _p; }
    virtual double V() const { return Model::extensiveProperty<&Model::v>(); }

    virtual Objective_t Y(size_t idx) const {
      switch(idx){
      case 0:
	return Objective_t::T;
      case 1:
	return Objective_t::p;
      default:
	stator_throw() << "Invalid variable index";
      };
    }

    void set(Objective_t val, double target, Objective_t constant) {
      if (val == Objective_t::T && constant == Objective_t::p)
	{ _T = target; return; }
      
      if (val == Objective_t::p && constant == Objective_t::T)
	{ _p = target; return; }
      
      stator_throw() << "Failed! set(" << val << ", " << target << ", " << constant << ")";
    }
    
    virtual double h(std::string molID) const {
      return chemPot(molID) + T() * s(molID);
    }
      
    virtual double u(std::string molID) const {
      return h(molID) - p() * v(molID);
    }

    virtual double a(std::string molID) const {
      return chemPot(molID) - p() * v(molID);
    }

    virtual double Y2() const { return _p; }
    virtual double& Y2() { return _p; }
    virtual double dfdT(Objective_t f) const {
      switch (f) {
      case Objective_t::G:
	return -S(); 
      case Objective_t::A:
	return -(S() + p() * Alpha() * V());
      case Objective_t::H:
	return Cp();
      case Objective_t::negS:
	return -Cp() / T();
      case Objective_t::S:
	return Cp() / T();
      case Objective_t::U:
	return Cp() - Alpha() * p() * V();
      case Objective_t::p:
	return 0;
      case Objective_t::V:
	return Alpha() * V();
      case Objective_t::T:
	return 1;
      default:
	stator_throw() << "Unhandled Objective_t " << f;
      }
    }
      
    virtual double dfdY2(Objective_t f) const {
      switch (f) {
      case Objective_t::G:
	return V();
      case Objective_t::H:
	return V() * (1 - T() * Alpha());
      case Objective_t::A:
	return Beta() * p() * V();
      case Objective_t::negS:
	return Alpha() * V();
      case Objective_t::S:
	return -Alpha() * V();
      case Objective_t::U:
	return V() * (Beta() * p() - Alpha() * T());
      case Objective_t::p:
	return 1;
      case Objective_t::V:
	return - Beta() * V();
      case Objective_t::T:
	return 0;
      default:
	stator_throw() << "Unhandled Objective_t " << f;
      }
    }
      
    virtual double dfdNi(Objective_t p, std::string molID) const {
      switch (p) {
      case Objective_t::G:
	return chemPot(molID);
      case Objective_t::H:
	return h(molID);
      case Objective_t::A:
	return a(molID);
      case Objective_t::negS:
	return -s(molID);
      case Objective_t::S:
	return s(molID);
      case Objective_t::U:
	return u(molID);
      case Objective_t::p:
	return 0;
      case Objective_t::V:
	return v(molID);
      case Objective_t::T:
	return 0;
      default:
	stator_throw() << "Unhandled Objective_t " << p;
	break;
      }
    }
      
  protected:
    double _p;
  };

  /*! \brief Ideal gas model using the variable set
    \f$\vec{X}=\left\{T,p,\left\{N_i\right\}^{N_c}\right\}\f$.
  */
  class ModelIdealGasTp : public ModelTp {
  public:
    ModelIdealGasTp(shared_ptr<Database> db, Components components, const double T, const double p):
      ModelTp(db, components, T, p, "idealgas")
    {}
      
    virtual double v(std::string molID) const {
      return _db->R() * T() / p();
    }

    virtual double s(std::string molID) const {
      //Get the temperature-only contribution
      const auto data = getGdata(molID, T());

      //Calculate the mixing entropy
      //Prevent divide by zero for empty phases
      double N = this->N();
      N += (N == 0);
      //Convert infinite mixing contributions to the maximum representation
      const double x = std::max((*this)[molID] / N, std::numeric_limits<double>::min());

      return -data.df - _db->R() * (std::log(x) + std::log(_p / data.p));
    }

    virtual double chemPot(std::string molID) const {
      //Get the pure temperature polynomial contribution
      const auto data = getGdata(molID, T());
      return data.f - T() * (s(molID) + data.df);
    }
    virtual double Cp() const {
      double sum(0);
      for (const auto& comp : *this)
	sum += comp.second * cp_ig0(comp.first);
      return sum;
    }

    virtual double Alpha() const {
      return 1.0 / T();
    }

    virtual double Beta() const {
      return 1.0 / _p;
    }

    /*! \brief Output a text string representation of this phase.
     */
    virtual std::string str() const {
      std::ostringstream os;
      os << "IdealGas{ T:" << T() << "K, p:" << p()/1e5 << "bar | " << Components::str();
      return os.str()+" }";
    }  
  };
    
  /*! \brief A base class for phases described using
    \f$\vec{X}=\left\{T,V,\left\{N_i\right\}^{N_c}\right\}\f$.
  */
  class ModelTV : public ModelT {
  public:
    ModelTV(shared_ptr<Database> db, Components components, const double T, const double V, const std::string type):
      ModelT(db, components, T, type),
      _V(V)
    {}

    virtual Objective_t Y(size_t idx) const {
      switch(idx){
      case 0:
	return Objective_t::T;
      case 1:
	return Objective_t::V;
      default:
	stator_throw() << "Invalid variable index";
      };
    }

    void set(Objective_t val, double target, Objective_t constant) {
      if (val == Objective_t::T && constant == Objective_t::V)
	{ _T = target; return; }
      
      if (val == Objective_t::V && constant == Objective_t::T)
	{ _V = target; return; }
      
      stator_throw() << "Failed! set(" << val << ", " << target << ", " << constant << ")";
    }
    
    //Functions which need to be overloaded by base classes
    /*! \brief Partial volar Helmholtz free energy. */
    virtual double volar_a(std::string molID) const = 0;

    /*! \brief Partial volar internal energy. */
    virtual double volar_u(std::string molID) const = 0;

    /*! \brief Isochoric heat capacity.*/
    virtual double CV() const = 0;

    /*! \brief \f$(\partial p/\partial V)_{T,\left\{N_{i}\right\}}\f$.*/
    virtual double dpdV() const = 0;

    /*! \brief \f$(\partial p/\partial T)_{V,\left\{N_{i}\right\}}\f$.*/
    virtual double dpdT() const = 0;

    /*! \brief \f$(\partial p/\partial N_i)_{T,V,\left\{N_{j\neq i}\right\}}\f$.*/
    virtual double dpdNi(std::string molID) const = 0;

    //This is still undefined from the base class
    //virtual double p() const = 0;

    virtual double v(std::string molID) const { return - dpdNi(molID) / dpdV();}
    virtual double V() const { return _V; }
    virtual double Cp() const { return CV() - T() * std::pow(dpdT(), 2) / dpdV(); }

    virtual double chemPot(std::string molID) const { return volar_a(molID); }
    virtual double a(std::string molID) const { return volar_a(molID) - p() * v(molID); }
    /*! \brief Partial volar Gibbs free energy. */
    virtual double volar_g(std::string molID) const { return chemPot(molID) - v(molID) * V() * dpdV(); }
    virtual double u(std::string molID) const { return volar_u(molID) + v(molID) * (T() * dpdT() - p()); }
    virtual double s(std::string molID) const { return (u(molID) - a(molID)) / T(); }
    /*! \brief Partial volar entropy. */
    virtual double volar_s(std::string molID) const { return s(molID) - v(molID) * dpdT(); }
    virtual double h(std::string molID) const { return u(molID) + p() * v(molID); }
    /*! \brief Partial volar enthalpy. */
    virtual double volar_h(std::string molID) const { return h(molID) - v(molID) * (V() * dpdV() + T() * dpdT()); }
    virtual double Alpha() const { return - dpdT() / (dpdV() * V()); }
    virtual double Beta() const { return - 1.0 / (dpdV() * V()); }

    virtual double dfdT(Objective_t f) const {
      switch (f) {
      case Objective_t::G:
	return V() * dpdT() - S();
      case Objective_t::A:
	return -S();
      case Objective_t::H:
	return CV() + V() * dpdT();
      case Objective_t::negS:
	return -CV() / T();
      case Objective_t::S:
	return CV() / T();
      case Objective_t::U:
	return CV();
      case Objective_t::V:
	return 0;
      case Objective_t::p:
	return dpdT();
      case Objective_t::T:
	return 1;
      default:
	stator_throw() << "Unhandled Objective_t " << f;
      }
      stator_throw() << "Should not reach here!";
    }
      
    virtual double dfdY2(Objective_t f) const {
      switch (f) {
      case Objective_t::G:
	return V() * dpdV();
      case Objective_t::A:
	return -p();
      case Objective_t::H:
	return V() * dpdV() + T() * dpdT();
      case Objective_t::negS:
	return -dpdT();
      case Objective_t::S:
	return dpdT();
      case Objective_t::U:
	return T() * dpdT() - p();
      case Objective_t::p:
	return dpdV();
      case Objective_t::V:
	return 1;
      case Objective_t::T:
	return 0;
      default:
	stator_throw() << "Unhandled Objective_t " << f;
      }
    }
      
    virtual double dfdNi(Objective_t p, std::string molID) const {
      switch (p) {
      case Objective_t::G:
	return volar_g(molID);
      case Objective_t::H:
	return volar_h(molID);
      case Objective_t::A:
	return volar_a(molID);
      case Objective_t::negS:
	return -volar_s(molID);
      case Objective_t::S:
	return volar_s(molID);
      case Objective_t::U:
	return volar_u(molID);
      case Objective_t::p:
	return dpdNi(molID);
      case Objective_t::V:
	return 0;
      case Objective_t::T:
	return 0;
      default:
	stator_throw() << "Unhandled Objective_t " << p;
      }
    }

  protected:
    double _V;
  };


  /*! \brief Ideal gas model using the variable set
    \f$\vec{X}=\left\{T,V,\left\{N_i\right\}^{N_c}\right\}\f$.
  */
  class ModelIdealGasTV : public ModelTV {
  public:
    ModelIdealGasTV(shared_ptr<Database> db, Components components, const double T, const double V):
      ModelTV(db, components, T, V, "idealgas")
    {}

    //Functions which need to be overloaded by base classes
    /*! \brief Partial volar Helmholtz free energy. */
    virtual double volar_a(std::string molID) const {
      //Get the pure temperature polynomial contribution
      const auto data = getGdata(molID, T());

      //Calculate the mixing entropy
      //Prevent divide by zero for empty phases
      double N = this->N();
      N += (N == 0);
      //Convert infinite mixing contributions to the maximum representation
      const double x = std::max((*this)[molID] / N, std::numeric_limits<double>::min());

      return data.f + _db->R() * T() * (std::log(N * _db->R() * T() / (V() * data.p)) + std::log(x));
    }

    /*! \brief Partial volar internal energy. */
    virtual double volar_u(std::string molID) const {
      const auto data = getGdata(molID, T());
      return data.f - _db->R() * T() - T() * data.df;
    }

    /*! \brief Isochoric heat capacity.*/
    virtual double CV() const {
      double sum_Cp=0;
      for (const auto& comp : *this)
	sum_Cp += comp.second * cp_ig0(comp.first);
	
      return sum_Cp - V() * dpdT();
    }

    /*! \brief \f$(\partial p/\partial V)_{T,\left\{N_{i}\right\}}\f$.*/
    virtual double dpdV() const { return - N() * _db->R() * T() / (V() * V()); }

    /*! \brief \f$(\partial p/\partial T)_{V,\left\{N_{i}\right\}}\f$.*/
    virtual double dpdT() const { return N() * _db->R() / V(); }
      
    /*! \brief \f$(\partial p/\partial N_i)_{T,V,\left\{N_{j\neq i}\right\}}\f$.*/
    virtual double dpdNi(std::string molID) const { return _db->R() * T() / V(); }

    virtual double p() const { return N() * _db->R() * T() / V(); }

    /*! \brief Output a text string representation of this phase.
     */
    virtual std::string str() const {
      std::ostringstream os;
      os << "IdealGas{ T:" << T() << ", V:" << V() << " | " << Components::str();
      return os.str()+" }";
    }
  };

  /*! \brief Incompressible phase model.
   */
  class ModelIncompressible : public ModelTp {
  public:
    ModelIncompressible(shared_ptr<Database> db, Components components, const double T, const double p, const std::string type, bool mixingEntropy=true):
      ModelTp(db, components, T, p, type),
      _mixingEntropy(mixingEntropy)
    {}
      
    virtual double v(std::string molID) const {
      return 0.0;
    }

    virtual double s(std::string molID) const {
      //Get the temperature-only contribution
      const auto data = getGdata(molID, T());

      double s = -data.df;
      if (_mixingEntropy) {
	//Calculate the mixing entropy
	//Prevent divide by zero for empty phases
	double N = this->N();
	N += (N == 0);
	//Convert infinite mixing contributions to the maximum representation
	const double x = std::max((*this)[molID] / N, std::numeric_limits<double>::min());
	s += - _db->R() * std::log(x);
      }
      return s;
    }

    virtual double chemPot(std::string molID) const {
      //Get the pure temperature polynomial contribution
      const auto data = getGdata(molID, T());

      const double no_mix_mu = data.f + v(molID) * (_p - data.p);

      if (!_mixingEntropy)
	return no_mix_mu;
      
      double N = this->N();
      N += (N == 0);
      //Convert infinite mixing contributions to the maximum representation
      const double x = std::max((*this)[molID] / N, std::numeric_limits<double>::min());
      
      return no_mix_mu + _db->R() * T() * std::log(x);
    }

    virtual double Cp() const {
      double sum(0);
      for (const auto& comp : *this)
	sum += comp.second * cp_ig0(comp.first);
      return sum;
    }

    virtual double Alpha() const {
      return 0;
    }

    virtual double Beta() const {
      return 0;
    }

    /*! \brief Output a text string representation of this phase.
     */
    virtual std::string str() const {
      std::ostringstream os;
      os << "ModelIncompressible{ T:" << T() << "K, p:" << p()/1e5 << "bar | " << Components::str();
      return os.str()+" }";
    }
      
  protected:
    bool _mixingEntropy;
  };

  class ModelIncompressibleSolid : public ModelIncompressible {
  public:
    ModelIncompressibleSolid(shared_ptr<Database> db, Components components, const double T, const double p):
      ModelIncompressible(db, components, T, p, "solid", false)
    {}
  };

  class ModelIncompressibleLiquid : public ModelIncompressible {
  public:
    ModelIncompressibleLiquid(shared_ptr<Database> db, Components components, const double T, const double p):
      ModelIncompressible(db, components, T, p, "liquid", true)
    {}
  };

  /*! \brief A container for a set of phases.
   */
  class MultiPhase : public std::vector<shared_ptr<Model> > {
  protected:
    typedef std::vector<shared_ptr<Model> > Base;
  public:
    using Base::Base;

    double T() const {
      if (empty())
	return 0;
      const double retT = front()->T();
      for (const auto& phase : *this)
	if (phase->T() != retT)
	  stator_throw() << "Multiphase has multiple temperatures, cannot determine which one you want.";
      return retT;
    }
    
    double p() const {
      double p = HUGE_VAL;
      for (const auto& phase : *this)
	p = std::min(p, phase->p());
      return (p != HUGE_VAL) ? p : 0;
    }
    
    double N() const { return sysProperty(&Model::N); }
    double M() const { return sysProperty(&Model::M); }
    double G() const { return sysProperty(&Model::G); }
    double H() const { return sysProperty(&Model::H); }
    double U() const { return sysProperty(&Model::U); }
    double V() const { return sysProperty(&Model::V); }
    double S() const { return sysProperty(&Model::S); }
    double A() const { return sysProperty(&Model::A); }
    double Cp() const { return sysProperty(&Model::Cp); }
    double Alpha() const { return sysProperty(&Model::Alpha); }
    double Beta() const { return sysProperty(&Model::Beta); }

    double f(Objective_t f) const {
      switch (f) {
      case Objective_t::G:
	return G();
      case Objective_t::H:
	return H();
      case Objective_t::negS:
	return -S();
      case Objective_t::S:
	return S();
      case Objective_t::U:
	return U();
      case Objective_t::T:
	return T();
      case Objective_t::A:
	return A();
      case Objective_t::p:
	return p();
      case Objective_t::V:
	return V();
      default:
	stator_throw() << "Unhandled Objective_t " << f;
      }
    }
    
  protected:
    template<class extProp>
    double sysProperty(extProp prop) const {
      double sum(0);
      for (const auto& p : *this)
	sum += ((*p).*prop)();
      return sum;
    }
  };

  /*! \brief A class which can calculate the equilibrium state of a
    set of phases.
  */
  class System : public MultiPhase {
  public:
    enum class Optimiser_t {SLSQP, NelderMead, IPopt};

    static Objective_t naturalPotential(Objective_t Y1, Objective_t Y2) {
	if ((Y1 == Objective_t::p) && (Y2 == Objective_t::T))
	  return Objective_t::G;
	else if ((Y1 == Objective_t::V) && (Y2 == Objective_t::T))
	  return Objective_t::A;
	else if ((Y1 == Objective_t::p) && (Y2 == Objective_t::S))
	  return Objective_t::H;
	else if ((Y1 == Objective_t::p) && (Y2 == Objective_t::H))
	  return Objective_t::negS;
	else if ((Y1 == Objective_t::V) && (Y2 == Objective_t::U))
	  return Objective_t::negS;
	else if ((Y1 == Objective_t::V) && (Y2 == Objective_t::S))
	  return Objective_t::U;
	else
	  stator_throw() << "Invalid combination of Y1 and Y2.";
    }

    void dfdX(Objective_t f, double* out, size_t N) const {
      if (f == Objective_t::p) {
	//The pressure derivative only concerns the first phase
	  size_t n_idx(0);
	  //First, the derivatives of the first phase's pressure by the molar amounts
	  for (const auto& mol : *front()) {
	    out[n_idx] = front()->dfdNi(Objective_t::p, mol.first);
	    ++n_idx;
	  }
	  out[n_idx] = front()->dfdY2(Objective_t::p);
	  ++n_idx;

	  //Skip every other phase's molar variables
	  for (; n_idx < N-1; ++n_idx)
	    out[n_idx] = 0;
	  
	  //The derivative of the system temperature
	  out[N-1] = front()->dfdT(Objective_t::p);
      } else if (f == Objective_t::T) {
	 //The temperature variables are all collapsed onto one variable.
	for (size_t i(0); i < N-1; ++i)
	  out[i] = 0;
	out[N-1] = 1; //dfdT
      } else {
	size_t i(0);
	double sum_dfdT(0);
	for (const auto& phase : *this) {
	  for (const auto& mol : *phase)
	    out[i++] = phase->dfdNi(f, mol.first);
	  out[i++] = phase->dfdY2(f);
	  sum_dfdT += phase->dfdT(f);
	}
	
	out[i++] = sum_dfdT; //dfdT

	if (i != N)
	  stator_throw() << "Mismatch in variable dimensions";	
      }
    }
    
    /* \brief Constructor for the System.

       \param max_eval The maximum number of Gibb's free energy
       evaluations to allow. Simulations are halted if the limit of
       precision is reached or this max iteration limit is
       reached. Generally, we hope that the limits of precision are
       reached during iterations, but this is a guard against
       problematic systems.
      
       \param constraintTol The tolerance to which the elemental or
       species constraints are upheld. This is a relative tolerance
       applied to each species/element. It should never be set too low,
       as the optimizer will reject perfectly good solution points if
       this constraint is violated!
      
       \param debug Print output on the progression of the
       optimisation.
      
       \param alg Select the algorithm used for optimisation.

       \param elemental Enable "reactive" simulations and allow any
       reordering of elements between species. If false, only transfer
       of species between phases is allowed.
    */
    System(Objective_t Y1 = Objective_t::p, Objective_t Y2 = Objective_t::T, bool reactive = false,
	   Optimiser_t alg = Optimiser_t::IPopt, double molarConstraintTol = 1e-8, size_t max_eval = 100,
	   bool debug=false, bool scale_problem = true, bool rank_reduce_material_constraints = true):
      _Y1(Y1), _Y2(Y2),
      _alg(alg), _reactive(reactive), _molarConstraintTol(molarConstraintTol),
      _max_eval(max_eval), _debug(debug), _scale_problem(scale_problem),
      _rank_reduce_material_constraints(rank_reduce_material_constraints)
    {}

    void setEnsemble(Objective_t Y1, Objective_t Y2) {
      _Y1 = Y1;
      _Y2 = Y2;
    }

    /*! \brief Calculate the equilibrium state of the system.
      
      \param Y1 Used to select if pressure P or volume V is held constant.
      \param Y2 Used to select if S, H, U, or T is held constant.
    */
    void equilibrate(double initialY1 = HUGE_VAL, double initialY2 = HUGE_VAL) {
      try {
	//Empty phases, or systems with zero components are counted as already at equilibrium
	if (this->empty() || (this->N() == 0))
	  return;

	if (dynamic_cast<const ModelExcess*>(this->front().get()))
	  stator_throw() << "ModelExcess cannot be the first phase in an equilibration";
	
	Eigen::IOFormat python(Eigen::FullPrecision, 0, ", ", ",\n", "[", "]", "[", "]");
	
	_potential = naturalPotential(_Y1, _Y2);
		
	//Make all temperatures the same as only one temperature
	//variable is inserted into the optimisation.
	//
	//We take the mass averaged temperature, using the mass as a
	//"thermal mass" (assuming everything has the same heat
	//capacity)
	
	double sumMT = 0;
	double sumM = 0;
	for (shared_ptr<Model>& phase : *this) {
	  if (dynamic_cast<const ModelExcess*>(phase.get())) continue;
	  sumMT += phase->M() * phase->T();
	  sumM += phase->M();
	}
	
	const double avgT = sumMT / sumM;
	
	for (shared_ptr<Model>& phase : *this)
	  if (!dynamic_cast<const ModelExcess*>(phase.get()))
	    phase->set(Objective_t::T, avgT, phase->Y(1));

	//Determine the dimensionality of the problem.
	//The system has a temperature variable
	_n = 1;

	//Each phase has either an independent p or V variable (unless
	//its an excess phase), and the component amounts
	for (const auto& p: *this)
	  _n += p->size() + (!dynamic_cast<const ModelExcess*>(p.get()));
	
	//Count the number of constraints.
	//Need to constrain either system p or V
	//Need to constrain either T, S, H, or U
	_m = 2;
	
	decltype(_molar_constraint_grad.jacobiSvd(Eigen::ComputeFullU)) svd;	
	_constraint_idx.clear();
	
	//Set up the elemental/molar constraints
	//First, identify what quantities to be constrained (each element or component)
	Components _constraints;
	for (const auto& phase : *this)
	  if (_reactive)
	    for (const auto& mol : *phase)
	      _constraints += mol.second * phase->db()->getComponent(mol.first).getElements();
	  else
	    _constraints += *phase;

	//Enumerate these constraints
	{
	  size_t m_idx(0);
	  for (const auto& c : _constraints)
	    _constraint_idx[c.first] = m_idx++;
	}

	//Now figure out the indexes of the other constraints
	const size_t mol_constraints = _constraint_idx.size();
	  
	_constraint_idx[std::string(_Y1)+"_{sys}"] = mol_constraints;
	_constraint_idx[std::string(_Y2)+"_{sys}"] = mol_constraints + 1;
	  
	//Now create the coefficient matrix for the set of linear
	//constraints
	_molar_constraint_grad = ConstraintGradMatrix::Zero(_constraints.size(), _n);
	size_t n_idx(0);
	if (_reactive)
	  for (const auto& phase : *this) {
	    for (const auto& mol : *phase) {
	      Components elements = phase->db()->getComponent(mol.first).getElements();
	      for (const auto& ele : elements) {
		size_t m_idx = _constraint_idx[ele.first];
		_molar_constraint_grad(m_idx, n_idx) += ele.second;
	      }
	      ++n_idx;
	    }

	    //Skip over the pressure/volume variable
	    ++n_idx;
	  }
	else
	  for (const auto& phase : *this) {
	    for (const auto& mol : *phase) {
	      size_t m_idx = _constraint_idx[mol.first];
	      _molar_constraint_grad(m_idx, n_idx) += 1;
	      ++n_idx;
	    }
	    //Skip over the pressure/volume variable
	    ++n_idx;
	  }

	if (_debug)
	  std::cout << "Molar constraint vector before elimination :" << std::endl << _molar_constraint_grad << std::endl;
	
	if (_rank_reduce_material_constraints) {
	  //Only a thin (rank-reduced) V is needed, but a full U
	  //required to restore to the original constraint vector
	  svd = _molar_constraint_grad.jacobiSvd(Eigen::ComputeThinV | Eigen::ComputeFullU);
	  
	  if (_debug) {
	    std::cout << "Singular value decomposition of constraint vec\n"
		      << "U=" << svd.matrixU().format(python) << std::endl
		      << "S=" << svd.singularValues().format(python) << std::endl
		      << "V*=" << svd.matrixV().transpose().format(python) << std::endl;
	  }
	  
	  _molar_constraint_grad = svd.matrixV().transpose();
	  _molar_constraint_grad.conservativeResize(svd.rank(), Eigen::NoChange);
	  if (svd.rank() != std::ptrdiff_t(_constraints.size())) {
	    if (_debug)
	      std::cout << "Identified " << _constraints.size() - svd.rank() << " redundant constraints" << std::endl;
	  }
	  
	  if (_debug)
	    std::cout << "Reduced constraints matrix :" << std::endl << _molar_constraint_grad << std::endl;
	} else if (_debug) {
	  std::cout << "Constraint elimination skipped" << std::endl;
	}
	
	//Add the elemental/molar constraints to the count
	_m += _molar_constraint_grad.rows();
	  	  
	if (_debug)
	  std::cout << "Building the g0 constraint initial values vector" << std::endl;
	buildConstraintGrad();

	//Construct the g0 vector
	_g0 = decltype(_g0)::Zero(_m);
	
	{
	  std::vector<double> g0_buf(_m, 0.0);
	  constraints(_m, g0_buf.data(), _n, nullptr, nullptr);
	  //Correct the Y1/Y2 entries
	  g0_buf[_m-2] = (initialY1 == HUGE_VAL) ? f(_Y1) : initialY1;
	  g0_buf[_m-1] = (initialY2 == HUGE_VAL) ? f(_Y2) : initialY2;
	  _g0 = MapMatrix(g0_buf.data(), _m, 1);
	  g0_buf.clear();
	  if (_debug) {
	    std::cout << "Checking constraints are zeroed" << std::endl;
	    constraints(_m, g0_buf.data(), _n, nullptr, nullptr);
	  }
	}
	
	switch (_alg) {
	case Optimiser_t::SLSQP:
	case Optimiser_t::NelderMead:
#ifdef SIMCEM_NLOPT
	  {
	    nlopt::opt optimizer;
	    if (_alg == Optimiser_t::SLSQP)
	      optimizer = nlopt::opt(nlopt::LD_SLSQP, _n);
	    else if (_alg == Optimiser_t::NelderMead) {
	      optimizer = nlopt::opt(nlopt::AUGLAG, _n);
	      nlopt::opt local_optimizer(nlopt::LN_NELDERMEAD, _n);
	      optimizer.set_local_optimizer(local_optimizer);
	    }
	    
	    if (_debug)
	      std::cout << "SETTING OBJECTIVE FUNCTION" << std::endl;
	    optimizer.set_min_objective(&NLoptObjective, static_cast<void*>(this));
	    optimizer.set_lower_bounds(0.0);
	    optimizer.set_maxeval(_max_eval);
	    if (_debug)
	      std::cout << "SETTING CONSTRAINT FUNCTION" << std::endl;
	    optimizer.add_equality_mconstraint(&NLoptConstraint, static_cast<void*>(this), std::vector<double>(_m, _molarConstraintTol));
	    double finalF;
	    std::vector<double> x(_n);
	    save_state(x.data(), _n);
	    try {
	      nlopt::result r = optimizer.optimize(x, finalF);
	      if (_debug) std::cout << "Optimiser halted successfully:" << NLOPT_exit_code(r) << std::endl;
	    } catch (nlopt::roundoff_limited) {
	      if (_debug)
		std::cout << "OPTIMISATION RESTARTING DUE TO ROUND-OFF LIMIT" << std::endl;
	    } catch (Ipopt::IpoptException& e) {
	      stator_throw() << "IPopt exception caught: " << e.Message();
	    } catch (std::exception& e) {
	      throw;
	    } catch (...) {
	      stator_throw() << "Unknown/unhandled exception";
	    }
	
	    restore_state(x.data(), _n);
	    buildConstraintGrad();
	    save_state(x.data(), _n);
	    //Restart the optimizer
	    try {
	      nlopt::result r = optimizer.optimize(x, finalF);	
	      if (_debug) std::cout << "Optimiser halted successfully:" << NLOPT_exit_code(r) << std::endl;
	    } catch (nlopt::roundoff_limited) {
	      if (_debug)
		std::cout << "RESTARTED OPTIMISATION ABORTED DUE TO ROUND-OFF LIMIT" << std::endl;
	    } catch (Ipopt::IpoptException& e) {
	      stator_throw() << "IPopt exception caught: " << e.Message();
	    } catch (std::exception& e) {
	      throw;
	    } catch (...) {
	      stator_throw() << "Unknown/unhandled exception";
	    }

	    if (_debug)
	      std::cout << "Optimised f = " << finalF << std::endl;
	
	    restore_state(x.data(), _n);
	    return;
	  }
#else
	  stator_throw() << "NLOPT support not available";	  
#endif
	case Optimiser_t::IPopt:
	  {
	    struct NLProblem : public Ipopt::TNLP {
	      NLProblem(System& sys, const Eigen::Matrix<double, Eigen::Dynamic, 1>& g0):
		_sys(sys),
		_g0(g0)
	      {
		x_init.resize(_sys._n);
		_sys.save_state(x_init.data(), _sys._n);
	      }
	      
	      System& _sys;
	      Eigen::Matrix<double, Eigen::Dynamic, 1> _g0;	      
	      std::vector<double> x_init;
	      std::vector<double> lagrangian_multipliers;
	      
	      virtual bool get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g, Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style) {
		n = _sys._n;
		m = _sys._m;
		nnz_jac_g = _sys._m * _sys._n; //Number of non-zero jacobian entries
		nnz_h_lag = _sys._n * (_sys._n-1) / 2 + _sys._n; //Hessian non-zeros (dense)
		index_style = Ipopt::TNLP::C_STYLE;
		return true;
	      }
	      
	      /** Method to return the bounds for my problem */
	      virtual bool get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u, Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u) {
		// the variables have lower bounds of 0
		// the variables have upper bounds of +inf 
		for (Ipopt::Index i = 0; i < n; i++) {
		  x_l[i] = 0.0;
		  x_u[i] = 2e19;
		}
		// Ipopt interprets any number greater than nlp_upper_bound_inf as
		// infinity. The default value of nlp_upper_bound_inf and nlp_lower_bound_inf
		// is 1e19 and can be changed through ipopt options.
		
		for (Ipopt::Index i = 0; i < m; i++)
		  g_l[i] = g_u[i] = _g0(i);
		
		return true;
	      }
	      
	      /** Method to return the starting point for the algorithm */
	      virtual bool get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x, bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U, Ipopt::Index m, bool init_lambda, Ipopt::Number* lambda) {
		// Here, we assume we only have starting values for x (not the dual variables)
		if (init_z || init_lambda || !init_x)
		  stator_throw() << "Invalid get_starting_point state";

		std::copy(x_init.data(), x_init.data() + n, x);
		return true;
	      }
	      
	      /** Method to return the objective value */
	      virtual bool eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number& obj_value) {
		try {
		  obj_value = _sys.objective(n, (new_x ? x : nullptr), nullptr);
		  return true;
		} catch (...) {
		  if (_sys._debug)
		    std::cout << "Evaluating f failed due to " <<_sys._error_details << std::endl;
		  return false;
		}
	      }
	      
	      /** Method to return the gradient of the objective */
	      virtual bool eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f) {
		try {
		  _sys.objective(n, (new_x ? x : nullptr), grad_f);
		  return true;
		} catch (...) {
		  if (_sys._debug)
		    std::cout << "Evaluating gradient of f failed due to " <<_sys._error_details << std::endl;
		  return false;
		}		  
	      }
	      
	      /** Method to return the constraint residuals */
	      virtual bool eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Number* g) {
		try {
		  _sys.constraints(m, g, n, (new_x ? x : nullptr), nullptr);
		  return true;
		} catch (...) {
		  if (_sys._debug)
		    std::cout << "Evaluating g failed due to " <<_sys._error_details << std::endl;
		  return false;
		}
	      }

	      /** Method to return:
	       *   1) The structure of the jacobian (if "values" is NULL)
	       *   2) The values of the jacobian (if "values" is not NULL)
	       */
	      virtual bool eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m,
				      Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index *jCol, Ipopt::Number* values) {
		try {
		  if ((iRow!=NULL) && (jCol!=NULL))
		    for (size_t i(0); i < _sys._m * _sys._n; ++i) {
		      //iRow is m_idx
		      //jCol is n_idx
		      iRow[i] = i / n;
		      jCol[i] = i % n;
		    }
		  
		  if (values != NULL)
		    _sys.constraints(m, nullptr, n, (new_x ? x : nullptr), values);
		  
		  return true;
		} catch (...) {
		  if (_sys._debug)
		    std::cout << "Evaluating jacobian of g failed due to " <<_sys._error_details << std::endl;
		  return false;
		}
	      }

	      /** Method to return:
	       *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
	       *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
	       */
	      virtual bool eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
				  Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* lambda,
				  bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow,
				  Ipopt::Index* jCol, Ipopt::Number* values) {
		return false;//No hessian available
	      }

	      /** @name Solution Methods */
	      //@{
	      /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
	      virtual void finalize_solution(Ipopt::SolverReturn status,
					     Ipopt::Index n, const Ipopt::Number* x, const Ipopt::Number* z_L, const Ipopt::Number* z_U,
					     Ipopt::Index m, const Ipopt::Number* g, const Ipopt::Number* lambda,
					     Ipopt::Number obj_value,
					     const Ipopt::IpoptData* ip_data,
					     Ipopt::IpoptCalculatedQuantities* ip_cq) {
		_sys.restore_state(x, n);

		lagrangian_multipliers = std::vector<double>(lambda, lambda + m);
	      }
	      //@}
	    };
	    
	    Ipopt::IpoptApplication app;
	    app.RethrowNonIpoptException(true);
	    app.Options()->SetNumericValue("tol", 1e-7);
	    app.Options()->SetIntegerValue("max_iter", _max_eval);
	    app.Options()->SetNumericValue("bound_relax_factor", 0);
	    app.Options()->SetStringValue("hessian_approximation", "limited-memory");
	    app.Options()->SetStringValue("expect_infeasible_problem", "yes");

	    if (_scale_problem)
	      app.Options()->SetStringValue("nlp_scaling_method", "gradient-based");
	    else
	      app.Options()->SetStringValue("nlp_scaling_method", "none");
	    
	    if (_debug) {
	      app.Options()->SetStringValue("derivative_test", "first-order");
	      app.Options()->SetIntegerValue("print_level", 5);
	      //_debug = false;
	    } else {
	      app.Options()->SetIntegerValue("print_level", 0);	      
	    }
	    Ipopt::ApplicationReturnStatus status = app.Initialize();
	    if (status != Ipopt::Solve_Succeeded)
	      stator_throw() << "*** Error during IpOpt initialization!";

	    Ipopt::SmartPtr<NLProblem> problem(new NLProblem(*this, _g0));
	    //IPopt likes to know what the constraints are set to, so
	    //it can estimate the problem scaling, so here we disable
	    //our own offset correction
	    _g0 = decltype(_g0)::Zero(_m);
	    
	    status = app.OptimizeTNLP(problem);
	    switch (status) {
	    case Ipopt::Solve_Succeeded:
	    case Ipopt::Solved_To_Acceptable_Level:
	    case Ipopt::Search_Direction_Becomes_Too_Small:
	    case Ipopt::Feasible_Point_Found: {
	      
	      //Map the lagrangians back into the non-rank-reduced space
	      Eigen::Matrix<double, Eigen::Dynamic, 1> lagrangians;
	      lagrangians.setZero(mol_constraints + 2);

	      if (_rank_reduce_material_constraints) {
		////Copy the reduced set of lagrangians out (svd.rank() <= mol_constraints)
		//for (int i(0); i < svd.rank(); ++i)
		//  lagrangians[i] = problem->lagrangian_multipliers[i];
		//
		////Copy the rest of the multipliers out
		//// BROKEN!!!
		////lagrangians[mol_constraints + 1] = problem->lagrangian_multipliers[svd.rank() + 1];
		////lagrangians[mol_constraints + 2] = problem->lagrangian_multipliers[svd.rank() + 2];
		////for (size_t i(mol_constraints); i < problem->lagrangian_multipliers.size(); ++i)
		////lagrangians[i] = problem->lagrangian_multipliers[i];
		//
		////Invert the non-zero singular values
		//Eigen::Matrix<double, Eigen::Dynamic, 1> Sinv;
		////Set all inverted singular values to 0, we do this as
		////there may be less than mol_constraints singular values
		////returned by the SVD calcluation (if the calc is "thin")
		//Sinv.setZero(mol_constraints);
		//for (int i(0); i < svd.rank(); ++i)
		//  if (svd.singularValues()(i) != 0)
		//    Sinv(i) = 1.0 / svd.singularValues()(i);
		//
		//if (_debug)
		//  std::cout << "lagrangians.head(mol_constraints) = " << lagrangians.head(mol_constraints).format(python) << std::endl
		//	    << "S^{-1} = " << Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(Eigen::DiagonalMatrix<double, Eigen::Dynamic>(Sinv)).format(python) << std::endl
		//	    << "U = " << svd.matrixU().format(python) << std::endl;
		//
		///////// BROKEN!!!
		//////// Working with Wahab I think the new lagrangians are this
		////lagrangians.head(mol_constraints) = lagrangians.head(mol_constraints) * Eigen::DiagonalMatrix<double, Eigen::Dynamic>(Sinv) * svd.matrixU().transpose();
		////Old one is here
		////lagrangians.head(mol_constraints) = svd.matrixU() * Eigen::DiagonalMatrix<double, Eigen::Dynamic>(Sinv) * lagrangians.head(mol_constraints);
	      } else {
		if (_debug)
		  std::cout << "Skipping lagrangian restoration as rank reduction is disabled" << std::endl;
		for (size_t i(0); i < problem->lagrangian_multipliers.size(); ++i)
		  lagrangians[i] = problem->lagrangian_multipliers[i];
	      }
	      
	      if (_debug)
		std::cout << "lagrangian array = " << lagrangians.format(python) << std::endl;
	      
	      for (const auto& p :  _constraint_idx)
		_lagrangians[p.first] = lagrangians[p.second];
	      
	      if (_debug)
		std::cout << "Found " << _constraint_idx.size() << " constraint indices" << std::endl
			  << "And " << _lagrangians.size() << " lagrangians" << std::endl;
	      break;
	    }
	    case Ipopt::Infeasible_Problem_Detected:
	      stator_throw() << "equilibrate() failed, IPopt: infeasible problem detected";
	    case Ipopt::Diverging_Iterates:
	      stator_throw() << "equilibrate() failed, IPopt: diverging iterates";
	    case Ipopt::User_Requested_Stop:
	      stator_throw() << "equilibrate() failed, IPopt: user requested stop";
	    case Ipopt::Maximum_Iterations_Exceeded:
	      stator_throw() << "equilibrate() failed, IPopt: maximum iterations exceeded";
	    case Ipopt::Restoration_Failed:
	      stator_throw() << "equilibrate() failed, IPopt: restoration failed";
	    case Ipopt::Error_In_Step_Computation:
	      stator_throw() << "equilibrate() failed, IPopt: error in step computation";
	    case Ipopt::Maximum_CpuTime_Exceeded:
	      stator_throw() << "equilibrate() failed, IPopt: CPU time exceeded";
	    case Ipopt::Not_Enough_Degrees_Of_Freedom:
	      stator_throw() << "equilibrate() failed, IPopt: not enough degrees of freedom";
	    case Ipopt::Invalid_Problem_Definition:
	      stator_throw() << "equilibrate() failed, IPopt: invalid problem definition";
	    case Ipopt::Invalid_Option:
	      stator_throw() << "equilibrate() failed, IPopt: invalid option";
	    case Ipopt::Invalid_Number_Detected:
	      stator_throw() << "equilibrate() failed, IPopt: invalid number detected";
	    case Ipopt::Unrecoverable_Exception:
	      stator_throw() << "equilibrate() failed, IPopt: unrecoverable exception";
	    case Ipopt::NonIpopt_Exception_Thrown:
	      stator_throw() << "equilibrate() failed, IPopt: Non-Ipopt exception thrown";
	    case Ipopt::Insufficient_Memory:
	      stator_throw() << "equilibrate() failed, IPopt: insufficient memory";
	    case Ipopt::Internal_Error:
	      stator_throw() << "equilibrate() failed, IPopt: internal IPopt error";
	    default:
	      stator_throw() << "*** Unknown Ipopt error " << status << " during IpOpt solution!";
	    };
	    return;
	  }
	};
      } catch (Ipopt::IpoptException& e) {
	stator_throw() << "IPopt exception caught: " << e.Message();
      } catch (std::exception& e) {
	throw;
      } catch (...) {
	stator_throw() << "Unknown/unhandled exception";
      }
    }

    std::string str() const {
      std::ostringstream os;
      os << "System { ";
      for (const auto& phase : *this)
	os << phase->str();
      os << " }";
      return os.str();
    }

    const std::map<std::string, double>& getLagrangians() const {
      return _lagrangians;
    }
    
  private:
    void save_state(double* const x, size_t N) const {
      size_t i = 0;
      for (const auto& phase : *this) {
	for (const auto& mol : *phase) {
	  x[i] = mol.second;
	  ++i;
	}
	if (!dynamic_cast<const ModelExcess*>(phase.get())) {
	  x[i] = phase->f(phase->Y(1));
	  ++i;
	}
      }

      x[i] = front()->T();
      ++i;
	
      if (i != N)
	stator_throw() << "Mismatch in variable dimensions";
    }

    void restore_state(const double* const x, size_t N) {      
      size_t i = 0;
      for (const auto& phase : *this) {
	for (auto& mol : *phase) {
	  mol.second = x[i];
	  ++i;
	}
	if (!dynamic_cast<const ModelExcess*>(phase.get())) {
	  phase->set(phase->Y(1), x[i], phase->Y(0));
	  ++i;
	}
      }
	
      const double newT = x[i];
      for (shared_ptr<Model>& phase : *this)
	if (!dynamic_cast<const ModelExcess*>(phase.get()))
	  phase->set(Objective_t::T, newT, phase->Y(1));
      ++i;

      if (i != N)
	stator_throw() << "Mismatch in variable dimensions";
    }

    void buildConstraintGrad() {
      //Build the molecular/elemental constraint gradients
      _constraint_grad = std::vector<double>(_m * _n, 0.0);
      for (ptrdiff_t m_idx(0); m_idx < _molar_constraint_grad.rows(); ++m_idx)
	for (size_t n_idx(0); n_idx < _n; ++n_idx)
	  _constraint_grad[m_idx * _n + n_idx] += _molar_constraint_grad(m_idx, n_idx);
    }

    double objective(const unsigned n, const double * const x, double * const grad) {
      try {
	if (x)
	  restore_state(x, n);

	const double f = MultiPhase::f(_potential);
	  
	if (grad)
	  dfdX(_potential, grad, n);

	if (_debug) {
	  std::cout << "Objective evaluation" << std::endl;
	  std::cout << "Objective = " << f << std::endl;
	  for (const auto& phase : *this)
	    std::cout << *phase << std::endl;

	  if (grad) {
	    std::cout << "Exact f differentials" << std::endl;
	    printObjectiveGradients(n, grad, _potential);
	    std::cout << "Numerical f differentials" << std::endl;
	    const auto dfdX_n = numerical_f_grad(_potential);
	    printObjectiveGradients(n, dfdX_n.data(), _potential);
	  }
	}
	
	return f;
      } catch (...) {
	throw; 
      }
   }
    
    void constraints(const unsigned m, double * const result, const unsigned n, const double * const x, double * const grad) {
      try {
	if (x)
	  restore_state(x, n);

	if (result) {
	  std::vector<double> x_state(n);
	  save_state(x_state.data(), n);
	  MapMatrix v_x(x_state.data(), n, 1);	  

	  //We only work with the molar constraints at this point
	  unsigned m_count = _molar_constraint_grad.rows();
	  {
	    //Here we only operate on the molar constraints
	    MapMatrix v_result(result, m_count, 1);
	    v_result = MapMatrix(_constraint_grad.data(), m_count, n) * v_x;
	  }

	  //Now we operate on the whole matrix
	  MapMatrix v_result(result, m, 1);
	  
	  //p or V constraint
	  switch(_Y1) {
	  case Objective_t::p:
	    result[m_count] = front()->p();
	    break;
	  case Objective_t::V:
	    result[m_count] = V();
	    break;
	  default:
	    stator_throw() << "Shouldn't reach here!";
	  }
	  ++m_count;

	  //T, S, H, or U constraint
	  result[m_count] = f(_Y2);
	  ++m_count;
	    
	  if (m_count != m)
	    stator_throw() << "Dimensions of problem have changed! m_current="<<m_count << " (old="<<m<<")";

	  v_result -= _g0;
	  
	  if (_debug) {
	    size_t m(0);
	    for (ptrdiff_t i(0); i < _molar_constraint_grad.rows(); ++i)
	      std::cout << "g_{" << i << "} = " << result[m++] << std::endl;
	    std::cout << "g_{"<<_Y1<<"} = " << result[m++] << std::endl;
	    std::cout << "g_{"<<_Y2<<"} = " << result[m++] << std::endl;
	  }
	}

	if (grad) {
	  //The molar constraints are linear, thus we just copy the
	  //pregenerated gradient matrix
	  size_t m_idx = _molar_constraint_grad.rows();
	  std::copy(_constraint_grad.data(), _constraint_grad.data() + m * n, grad);

	  //The p or V constraint
	  dfdX(_Y1, grad + m_idx * n, n);
	  ++m_idx;
	  
	  //The S, H, U, or T constraint
	  dfdX(_Y2, grad + m_idx * n, n);
	  ++m_idx;
	    
	  if (m_idx != m)
	    stator_throw() << "Dimensions of problem have changed! m_current="<<m_idx << " (old="<<m<<")";

	  if (_debug) {
	    std::cout << "Exact g differentials" << std::endl;
	    printConstraintGradients(m, n, grad);
	    std::cout << "Numerical g differentials" << std::endl;
	    const auto dgdX = numerical_constraint_grad();
	    printConstraintGradients(m, n, dgdX.data());
	  }
	}
#ifdef SIMCEM_NLOPT
      } catch (std::exception & e) {

	_error_details = e.what();
	throw nlopt::forced_stop();
#endif
      } catch (...) {
	stator_throw() << "Unknown/unhandled exception";
      }
    }
    
    void printStateVector(size_t n, const double* x) const {
      std::ios::fmtflags old_flags(std::cout.flags());
      std::cout << std::fixed << std::left << std::setprecision(8);
      size_t phase_idx(0);
      size_t idx(0);
      for (const auto& phase : *this) {
	for (const auto& mol : *phase)
	  std::cout << std::setw(10) << (mol.first.substr(0, 8))
		    << std::setw(10) << x[idx++] << std::endl;
	if (dynamic_pointer_cast<ModelTp>(phase))
	  std::cout << std::setw(10) << ("p("+std::to_string(phase_idx)+")")
		    << std::setw(10) << x[idx++] << std::endl;
	else if (dynamic_pointer_cast<ModelTV>(phase))
	  std::cout << std::setw(10) << ("V("+std::to_string(phase_idx)+")")
		    << std::setw(10) << x[idx++] << std::endl;
	else
	  std::cout << std::setw(10) << ("Y2("+std::to_string(phase_idx)+")")
		    << std::setw(10) << x[idx++] << std::endl;
	++phase_idx;
      }
      std::cout << std::setw(10) << "T_{sys}" << std::setw(10) << x[idx++] << std::endl;
      std::cout.flags(old_flags);
    }

    void printObjectiveGradients(size_t n, const double* grad, const Objective_t f) const {
      std::ios::fmtflags old_flags(std::cout.flags());
      std::cout << std::fixed << std::setprecision(4)
		<< std::setw(10) << std::left << " ";

      size_t phase_idx(0);
      for (const auto& phase : *this) {
	for (const auto& mol : *phase)
	  std::cout << std::setw(10) << ("d/d"+mol.first+"("+std::to_string(phase_idx)+")");
	if (dynamic_pointer_cast<ModelTp>(phase))
	  std::cout << std::setw(10) << ("d/dp("+std::to_string(phase_idx)+")");
	else if (dynamic_pointer_cast<ModelTV>(phase))
	  std::cout << std::setw(10) << ("d/dV("+std::to_string(phase_idx)+")");
	else
	  std::cout << std::setw(10) << ("d/dY_2("+std::to_string(phase_idx)+")");
	++phase_idx;
      }

      std::cout << std::setw(10) << "d/dT_{sys}" << std::endl
		<< std::setw(10) << f;

      size_t idx(0);
      for (const auto& phase : *this) {
	for (size_t i(0); i < phase->size(); ++i)
	  std::cout << std::setw(10) << grad[idx++];
	std::cout << std::setw(10) << grad[idx++];
      }

      std::cout << std::setw(10) << grad[idx++] << std::endl;

      std::cout.flags(old_flags);
      if (idx != n)
	stator_throw() << "Mismatch in variable dimensions";
    }
    
    void printConstraintGradients(size_t m, size_t n, const double* grad) const {
      std::ios::fmtflags f(std::cout.flags());
      std::cout << std::left;
      std::cout << std::setw(10) << "";
      //Print table header
      size_t phase_counter(0);
      for (const auto& phase : *this) {
	for (const auto& mol : *phase)
	  std::cout << std::setw(10) << (mol.first + "(" + std::to_string(phase_counter) + ")");

	if (dynamic_pointer_cast<ModelTp>(phase))
	  std::cout << std::setw(10) << ("p("+std::to_string(phase_counter)+")");
	else if (dynamic_pointer_cast<ModelTV>(phase))
	  std::cout << std::setw(10) << ("V("+std::to_string(phase_counter)+")");
	else
	  std::cout << std::setw(10) << ("Y_2("+std::to_string(phase_counter)+")");

	++phase_counter;
      }
      std::cout << std::setw(10) << "T_{sys}" << std::endl;

      //Now print the elemental/species constraints
      ptrdiff_t m_idx(0);
      for (; m_idx < _molar_constraint_grad.rows(); ++m_idx) {
	std::cout << std::setw(10) << " ";
	for (size_t j(0); j < n; ++j)
	  std::cout << std::setw(10) << grad[m_idx * n + j];
	std::cout << std::endl;
      }

      std::cout << std::setw(10) << _Y1;
      for (size_t j(0); j < n; ++j)
	std::cout << std::setw(10) << grad[m_idx * n + j];
      std::cout << std::endl;
      ++m_idx;
      
      std::cout << std::setw(10) << _Y2;
      for (size_t j(0); j < n; ++j)
	std::cout << std::setw(10) << grad[m_idx * n + j];
      std::cout << std::endl;
      ++m_idx;
      std::cout.flags(f);
    }

    std::vector<double> numerical_f_grad(const Objective_t potential, const double diff = 1.0e-6) {
      std::vector<double> dfdX(_n);
      std::vector<double> x(_n);
      save_state(x.data(), _n);
      const double f_orig = MultiPhase::f(potential);

      for (size_t n_idx(0); n_idx < _n; ++n_idx) {
	const double oldx = x[n_idx];
	double dX = diff * (x[n_idx] + (x[n_idx]==0));
	x[n_idx] += dX;
	restore_state(x.data(), _n);

	const double f_new = MultiPhase::f(potential);
	dfdX[n_idx] = (f_new - f_orig) / dX;
	x[n_idx] = oldx;
      }

      restore_state(x.data(), _n);
      return dfdX;
    }

    std::vector<double> numerical_constraint_grad(const double diff = 1.0e-6) {
      std::vector<double> dgdX(_n * _m);

      //Save original state
      std::vector<double> x(_n);
      save_state(x.data(), _n);
      bool debug_flag = false; //Have to disable debug!
      std::swap(_debug, debug_flag);

      //Calculate original values
      std::vector<double> g_orig(_m);
      constraints(_m, g_orig.data(), _n, nullptr, nullptr);
      std::vector<double> g_new(_m);
      
      for (size_t n_idx(0); n_idx < _n; ++n_idx) {
	const double oldx = x[n_idx];
	double dX = diff * x[n_idx];
	if (dX == 0) dX = 1e-12;
	
	x[n_idx] += dX;
	restore_state(x.data(), _n);

	constraints(_m, g_new.data(), _n, nullptr, nullptr);

	for (size_t m_idx(0); m_idx < _m; ++m_idx)
	  dgdX[m_idx * _n + n_idx] = (g_new[m_idx] - g_orig[m_idx]) / dX;
	x[n_idx] = oldx;
      }

      //Restore original state
      std::swap(_debug, debug_flag);
      restore_state(x.data(), _n);
      return dgdX;
    }

#ifdef SIMCEM_NLOPT
    static std::string NLOPT_exit_code(nlopt::result v) {
      switch (v) {
      case nlopt::SUCCESS: return "Success";
      case nlopt::STOPVAL_REACHED: return "Stopval reached";
      case nlopt::FTOL_REACHED: return "F tolerance reached";
      case nlopt::XTOL_REACHED: return "X tolerance reached";
      case nlopt::MAXEVAL_REACHED: return "Max evaluations reached";
      case nlopt::MAXTIME_REACHED: return "Max time reached";
      default:
	return "Unknown exit code! ";
      }
    }
#endif
    
    static double NLoptObjective(unsigned n, const double *x, double *grad, void *system) {
      return reinterpret_cast<simcem::System*>(system)->objective(n, x, grad);
    }

    static void NLoptConstraint(unsigned m, double *result, unsigned n, const double *x, double *grad, void *system) {
      reinterpret_cast<simcem::System*>(system)->constraints(m, result, n, x, grad);
    }
    
    Objective_t _Y1;
    Objective_t _Y2;
    Optimiser_t _alg;
    Objective_t _potential = Objective_t::G;
    bool _reactive;
    std::vector<double> _constraint_grad;
    std::string _error_details;
      
    size_t _n;
    size_t _m;

    double _molarConstraintTol;
    size_t _max_eval;
    bool _debug;
    bool _scale_problem;
    bool _rank_reduce_material_constraints;
    
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> ConstraintGradMatrix;
    ConstraintGradMatrix _molar_constraint_grad;

    typedef Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > MapMatrix;
    
    Eigen::Matrix<double, Eigen::Dynamic, 1> _g0;
    std::map<std::string, size_t> _constraint_idx;
    std::map<std::string, double> _lagrangians;
  };

  inline std::ostream& operator<<(std::ostream& os, const System& p) {
    return os << p.str();
  }

  namespace trans {

    /*! \brief NASA CEA calculation of the viscosity and thermal
        conductivity of the gases.

      As outlined in the reference NASA_CEA_Analysis.
     */
    inline
    std::pair<double, double> NASA_transport(const shared_ptr<Model> phase) {
      //Get a cache of the ID, moles, mass, viscosity, and thermal conductivity from each isobar.
      std::vector<std::tuple<std::string, double, double, double, double> > dataCache;
      dataCache.reserve(phase->size());

      const Database& db = *(phase->db());

      double Ntotal = 0;
      for (const Components::value_type& i : *phase) {
	IsobarRange isobars(db, i.first,
			    {[](const Phase& phase) { return phase._type != "idealgas"; }},
			    {[](const shared_ptr<Isobar>& isobar) {
				return
				  ((isobar->getCurve()->getVariable() != Objective_t::dynvisc)
				   && (isobar->getCurve()->getVariable() != Objective_t::thermcond))
				  || (isobar->getCurve()->_source != "NASA_GLENN"); }}
			    );
	double visc = 0;
	double thermcond = 0;
	sym::Var<sym::vidx<'T'>> varT;
	for (const shared_ptr<Isobar>& isobar: isobars)
	  {
	    const auto& curve = isobar->getCurve();
	    if (curve->inRange(phase->T())) {
	      Ntotal += i.second;
	      if (isobar->getCurve()->getVariable() == Objective_t::dynvisc)
		visc = sym::fast_sub(isobar->getCurve()->getFunction(), varT = phase->T());
	      if (isobar->getCurve()->getVariable() == Objective_t::thermcond)
		thermcond = sym::fast_sub(isobar->getCurve()->getFunction(), varT = phase->T());
	    }
	  }

	if (visc && thermcond)
	  dataCache.push_back(decltype(dataCache)::value_type(i.first, i.second, db.getComponent(i.first).mass(), visc, thermcond));
      }

      if (Ntotal <= 0.9 * phase->N())
	stator_throw() << "Only " << Ntotal / (phase->N() + (phase->N()==0))*100 << " mol% of " << *phase << " has viscosity data!";
      
      double visMix = 0, thermMix = 0;
      for (size_t i(0); i < dataCache.size(); ++i) {
	const double& N_i = std::get<1>(dataCache[i]);
	const double& m_i = std::get<2>(dataCache[i]);
	const double& nu_i = std::get<3>(dataCache[i]);
	const double& therm_i = std::get<4>(dataCache[i]);
	double visc_denom = N_i;
	double therm_denom = N_i;
	for (size_t j(0); j < dataCache.size(); ++j) {
	  if (i == j) continue;
	  const double& N_j = std::get<1>(dataCache[j]);
	  const double& m_j = std::get<2>(dataCache[j]);
	  const double& nu_j = std::get<3>(dataCache[j]);

	  const double Phi_ij = 0.25 * std::pow(1 + std::sqrt((nu_i / nu_j) * std::sqrt(m_j / m_i)), 2) * std::sqrt( 2 * m_j / (m_i + m_j));

	  const double Psi_ij = Phi_ij * (1 + 2.41 * (m_i - m_j) * (m_i - 0.142 * m_j) / std::pow(m_i + m_j, 2));
	  
	  visc_denom += N_j * Phi_ij;
	  therm_denom += N_j * Psi_ij;
	}

	visMix += N_i * nu_i / visc_denom;
	thermMix += N_i * therm_i / therm_denom;
      }

      return std::make_pair(visMix, thermMix);
    }
  }
  
  class SystemNASACEA : public System {
    SystemNASACEA(shared_ptr<Database> db, double T, double p, bool reactive = false) : System(reactive) {
      //Create the gas and solid phases
      push_back(shared_ptr<Model>(new ModelIdealGasTp(db, Components(), T, p)));
      push_back(shared_ptr<Model>(new ModelIncompressible(db, Components(), T, p, "solid")));
    }
  };
}
