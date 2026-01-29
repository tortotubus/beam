{
  description = "ELFF (source-tree)";

  # Optional: helps CI/users if they haven't enabled flakes globally
  nixConfig = {
    experimental-features = [ "nix-command" "flakes" ];
  };

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
  inputs.tortotubus.url = "github:tortotubus/nur-packages";
  inputs.tortotubus.inputs.nixpkgs.follows = "nixpkgs";

  outputs = { self, nixpkgs, tortotubus}:
    let
      systems = [ "x86_64-linux" "aarch64-linux" ];
      forAllSystems = f: nixpkgs.lib.genAttrs systems (system: f system);
    in
    {
      overlays.default = final: prev: {
        elff = final.callPackage ./elff.nix { 
          pacific = tortotubus.packages.${final.system}.pacific;
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
            ];
          };
        });

      # Nice: `nix flake check` builds it
      checks = forAllSystems (system: {
        elff = self.packages.${system}.elff;
      });
    };
}
