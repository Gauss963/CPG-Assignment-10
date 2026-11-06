#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
build_root="$RUNNER_TEMP/assignment-build"
object_dir="$build_root/objects"
release_dir="$RUNNER_TEMP/release-assets"
bundle_dir="$build_root/macos-bundle"
pgplot_source="$build_root/pgplot-source"
pgplot_build="$build_root/pgplot-build"
flang="$(brew --prefix flang)/bin/flang"
pgplot_revision="007707a71b18e6722be1ee1551dc3f654d1be511"
mkdir -p "$object_dir" "$release_dir" "$bundle_dir/bin" "$bundle_dir/pgplot"

source_count=0
while IFS= read -r -d '' source_file; do
  source_count=$((source_count + 1))
  "$flang" -c -J "$object_dir" -I "$object_dir" \
    -I "$(dirname "$source_file")" "$source_file" -o "$object_dir/$source_count.o"
done < <(
  find "$repo_root" -type f \
    \( -iname '*.f' -o -iname '*.for' -o -iname '*.f77' -o \
       -iname '*.f90' -o -iname '*.f95' -o -iname '*.f03' -o \
       -iname '*.f08' \) \
    -not -path '*/.git/*' -print0
)
if [ "$source_count" -eq 0 ]; then
  echo "No Fortran source files were found." >&2
  exit 1
fi

archive="$build_root/pgplot.tar.gz"
curl --fail --location --retry 3 \
  "https://github.com/astroumd/pgplot/archive/$pgplot_revision.tar.gz" -o "$archive"
mkdir -p "$pgplot_source" "$pgplot_build"
tar -xzf "$archive" --strip-components=1 -C "$pgplot_source"
cp "$pgplot_source/sys_macosx/gfortran_gcc_64.conf" "$pgplot_source/sys_macosx/flang_clang.conf"
sed -i '' \
  -e 's/FCOMPL="gfortran"/FCOMPL="flang"/' \
  -e 's/FFLAGC="[^"]*"/FFLAGC="-O2 -fPIC"/' \
  -e 's/FFLAGD="[^"]*"/FFLAGD="-O2"/' \
  -e 's/CCOMPL="gcc"/CCOMPL="clang"/' \
  -e 's/CFLAGC="[^"]*"/CFLAGC="-Wall -fPIC -DPG_PPU -O2"/' \
  -e 's/CFLAGD="[^"]*"/CFLAGD="-Wall -O2"/' \
  -e 's|LIBS="[^"]*"|LIBS=""|' \
  "$pgplot_source/sys_macosx/flang_clang.conf"
cp "$pgplot_source/drivers.list" "$pgplot_build/drivers.list"
sed -i '' -E 's/^! (PSDRIV [1-4] )/  \1/' "$pgplot_build/drivers.list"
(
  cd "$pgplot_build"
  "$pgplot_source/makemake" "$pgplot_source" macosx flang_clang
  make -j"$(sysctl -n hw.ncpu)" lib grfont.dat
)
cp "$pgplot_build/grfont.dat" "$pgplot_source/rgb.txt" \
  "$pgplot_source/copyright.notice" "$bundle_dir/pgplot/"

while IFS='|' read -r target sources link_flags; do
  [[ -z "$target" || "$target" == \#* ]] && continue
  read -r -a source_files <<< "$sources"
  binary="$release_dir/$target-macos-arm64"
  (
    cd "$repo_root"
    if [[ -n "${link_flags:-}" ]]; then
      read -r -a extra_flags <<< "$link_flags"
      "$flang" -J "$object_dir" -I "$object_dir" "${source_files[@]}" \
        "$pgplot_build/libpgplot.a" "${extra_flags[@]}" -o "$binary"
    else
      "$flang" -J "$object_dir" -I "$object_dir" "${source_files[@]}" \
        "$pgplot_build/libpgplot.a" -o "$binary"
    fi
  )
  cp "$binary" "$bundle_dir/bin/$target"
done < "$repo_root/.github/scripts/release-targets.txt"

cat > "$bundle_dir/run.sh" <<'RUN'
#!/bin/sh
set -eu
if [ "$#" -lt 1 ]; then
  echo "Usage: $0 PROGRAM [arguments ...]" >&2
  exit 2
fi
bundle_dir="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
program="$1"
shift
export PGPLOT_DIR="$bundle_dir/pgplot"
exec "$bundle_dir/bin/$program" "$@"
RUN
chmod +x "$bundle_dir/run.sh"
cat > "$bundle_dir/README.txt" <<'README'
macOS arm64 Fortran executables built with LLVM Flang.

Run a program from the repository directory containing its required input data:
  /path/to/macos-bundle/run.sh PROGRAM_NAME [arguments ...]

The launcher sets PGPLOT_DIR for the bundled PostScript-capable PGPLOT runtime.
PGPLOT copyright and redistribution terms are included in pgplot/copyright.notice.
README
ditto -c -k --sequesterRsrc --keepParent "$bundle_dir" \
  "$release_dir/$(basename "$repo_root")-macos-arm64.zip"
