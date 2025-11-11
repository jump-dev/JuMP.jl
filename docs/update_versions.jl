#  Copyright 2023, Oscar Dowson and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at https://mozilla.org/MPL/2.0/.

import Pkg
import TOML

"""
    list_changes()

Print a list of discrepancies between the Julia registry and the `packages.toml`
and `Project.toml` files.
"""
function list_changes()
    data = Pkg.Registry.uncompress_registry(
        expanduser("~/.julia/registries/General.tar.gz"),
    )
    println("docs/packages.toml")
    for (name, pkg_data) in TOML.parsefile("docs/packages.toml")
        name = get(pkg_data, "package", name)
        versions_file = data["$(uppercase(first(name)))/$(name)/Versions.toml"]
        versions = TOML.parse(versions_file)
        latest_version = maximum(VersionNumber.(keys(versions)))
        latest_version_s = "v$latest_version"
        if pkg_data["rev"] != latest_version_s
            println("$name")
            println("  packages.toml: ", pkg_data["rev"])
            println("       registry: ", latest_version_s)
            println()
        end
    end
    println("docs/Project.toml")
    for (name, version) in TOML.parsefile("docs/Project.toml")["compat"]
        versions_file = data["$(uppercase(first(name)))/$(name)/Versions.toml"]
        versions = TOML.parse(versions_file)
        latest_version = maximum(VersionNumber.(keys(versions)))
        latest_version_s = "$latest_version"
        if replace(version, "=" => "") != latest_version_s
            println("$name")
            println("  Project.toml: ", version)
            println("      registry: ", latest_version_s)
            println()
        end
    end
    return
end

list_changes()
