using(DelimitedFiles)

"""
This script takes a origin directory path, file name in the origin directory, 
a destiny directory path, and file name in the destiny directory. 
If the filename in the origin directory exists in the destiny
directory, then the two files are merged in the destiny directory.
If the filename in the origin directory does not exist in the destiny 
directory, then a file with that name is created in the destiny directory.
"""
function merge_files(origin_dir_name::String, 
    origin_file_name::String,
    destiny_dir_name::String,
    destiny_file_name::String)
    #check if origin_file_name exists in destiny directory   
    exists = isfile("..data/"*destiny_dir_name*"/"*origin_file_name)
    if exists
        old_file = readdlm("../data/"*destiny_dir_name*"/"*destiny_file_name) #big file
        new_file = readdlm("../data/"*origin_dir_name*"/"*origin_file_name,
                           skipstart=1) #file to merge
        #merge with the existing file
        big_file = vcat(old_file, new_file)
        #save the bigger file
        writedlm("../data/"*destiny_dir_name*"/"*destiny_file_name, big_file)
    else
        #load new file
        new_file = readdlm("../data/"*origin_dir_name*"/"*origin_file_name)
        #create the file and write to it
        io = open("../data/"*destiny_dir_name*"/"*destiny_file_name, "w")
        writedlm("../data/"*destiny_dir_name*"/"*destiny_file_name, new_file)
    end                 
end