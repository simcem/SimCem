{ pkgs ? import <nixpkgs> {} }:
{
 simcem = pkgs.callPackage ./derivation.nix {};
}
