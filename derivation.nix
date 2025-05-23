# see https://unix.stackexchange.com/questions/717168/how-to-package-my-software-in-nix-or-write-my-own-package-derivation-for-nixpkgs
# To install `nix-env -u -f default.nix`
# To develop `nix-shell` (will build the shell with dependencies)
# To test build `nix-build`
{ pkgs, python3, ipopt, eigen, boost }:
python3.pkgs.buildPythonPackage rec {
  name = "simcem";
  src = ./.;
  pyproject = true;

  dontUseCmakeConfigure = true;

  # Build-time dependencies for the package
  nativeBuildInputs = with pkgs; [
    pkg-config
    cmake
    ninja
  ];
    
  # Runtime dependencies for the package
  buildInputs = with python3.pkgs; [
    python
    scikit-build-core
    setuptools
    cmake
  ] ++ [
    ipopt
    eigen
    boost
  ];
  
  dependencies = with python3.pkgs; [
    numpy
    scipy
    pytest
    pkgconfig
  ];
}
