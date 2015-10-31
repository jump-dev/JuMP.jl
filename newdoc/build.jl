import Mustache

all_md_files = filter(f->contains(f,".md"), readdir())

manpage_temp = readall(joinpath("templates","manpage.html"))

for md_file in all_md_files
    md_content   = readall(md_file)
    html_content = Markdown.html(Markdown.parse(md_content))
    out_content  = Mustache.render(manpage_temp,
                        Dict("CONTENT" => html_content))
    root_name, _ = splitext(md_file)
    output_name  = joinpath("build", string(root_name, ".html"))
    open(output_name, "w") do fp
        print(fp, out_content)
    end
end
