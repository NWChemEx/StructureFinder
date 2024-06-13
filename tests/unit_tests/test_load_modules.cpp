#include <catch2/catch.hpp>
#include <structurefinder/structurefinder.hpp>

TEST_CASE("load_modules") {
    pluginplay::ModuleManager mm;
    structurefinder::load_modules(mm);
}
