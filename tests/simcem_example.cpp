#include <simcem/simcem.hpp>

int main() {
    const auto db = simcem::Database::create_from_file("../pysrc/simcem/free_database.xml");

    // Perform a combustion calculation using the database
    const auto gas = std::shared_ptr<simcem::ModelIdealGasTp>(new simcem::ModelIdealGasTp(db, {{"CH4",1}, {"O2",2}, {"N2",7.52}}, 298.15, 1e5));
    auto system = simcem::System(simcem::Objective_t::p, simcem::Objective_t::T, false,
                                simcem::System::Optimiser_t::IPopt, 1e-8, 100, true, true, true);
    system.push_back(gas);
    system.equilibrate();

    //Output results
    for (const auto& phase : system) {
        std::cout << phase->str() << std::endl;
    }
    return 0;
}