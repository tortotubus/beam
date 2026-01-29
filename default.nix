{ pkgs ? import <nixpkgs> {} }:
let
  tortotubusSrc = builtins.fetchTarball "https://github.com/tortotubus/nur-packages/archive/refs/heads/main.tar.gz";
  tortotubusOverlay = import (tortotubusSrc + "/overlay.nix");
  pkgs' = import <nixpkgs> { overlays = [ tortotubusOverlay ]; };
in
pkgs'.callPackage ./elff.nix { }
