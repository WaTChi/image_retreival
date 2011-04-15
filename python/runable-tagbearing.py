# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="zhangz"
__date__ ="$Apr 11, 2011 4:16:34 PM$"

if __name__ == "__main__":
    print "TAGGING"

    import tags
    import util

    input="/media/DATAPART2/jz/tags.csv"
    intermediate_file_for_ge="/media/DATAPART2/jz/gepathTags.txt"
    output="/media/DATAPART2/jz/earthminetags2.csv"

    t = tags.TagCollection(input)
    t.output_to_gepath_format(intermediate_file_for_ge)
    util.generate_vector_tags(intermediate_file_for_ge, output)