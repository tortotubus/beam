{
  description = "ELFF (source-tree)";

  # Optional: helps CI/users if they haven't enabled flakes globally
  nixConfig = {
    experimental-features = [ "nix-command" "flakes" ];
  };

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
  inputs.basilisk = {
    url = "github:tortotubus/pacific-basilisk/17bf4c16f91abb6bfb0013501cabb80a645720df";
    flake = false;
  };

  outputs = { self, nixpkgs, basilisk }:
    let
      systems = [ "x86_64-linux" "aarch64-linux" ];
      forAllSystems = f: nixpkgs.lib.genAttrs systems (system: f system);
    in
    {
      overlays.default = final: prev: {
        elff = final.callPackage ./elff.nix {
          basiliskSrc = basilisk;
        };
      };

      packages = forAllSystems (system:
        let
          pkgs = import nixpkgs {
            inherit system;
            overlays = [ self.overlays.default ];
          };
        in
        {
          elff = pkgs.elff;
          default = pkgs.elff;
        });

      devShells = forAllSystems (system:
        let
          pkgs = import nixpkgs {
            inherit system;
            overlays = [ self.overlays.default ];
          };
        in
        {
          default = pkgs.mkShell {
            inputsFrom = [ pkgs.elff ];
            packages = with pkgs; [
              gdb
              strace
              ninja
              doxygen
              apptainer
              criterion
              gnuplot_qt
            ];
            shellHook = ''
              export GCOV="gcov -i"
            '';
          };
        });

      # Nice: `nix flake check` builds it
      checks = forAllSystems (system: {
        elff = self.packages.${system}.elff;
      });
    };
}
