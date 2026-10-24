{ pkgs ? import <nixpkgs> {} }:
let
  basiliskSrc = builtins.fetchTarball
    "https://github.com/tortotubus/pacific-basilisk/archive/17bf4c16f91abb6bfb0013501cabb80a645720df.tar.gz";
in
pkgs.callPackage ./elff.nix { inherit basiliskSrc; }
