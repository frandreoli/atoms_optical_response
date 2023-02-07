#
#This code calculates in real time the amount of RAM and CPU used by a target process
#It is intended to work on a Linux environment
#Tested on CentOS Linux 7 (Core)
#
const time_start = time()
#@time using HDF5
time_initialization = @timed using Plots
println("#Initialization time: ", time_initialization[2]," seconds")
ENV["GKSwstype"]="nul"

const file_name   = ARGS[1]
const process_id  = ARGS[2]
const max_hours   = parse(Float64,ARGS[3])
const cmd_ram = `ps -o rss $process_id`  #const cmd_ram = `pmap $process_id`#const cmd_ram=`ps -o vsz $process_id`
const cmd_cpu = `top -b -n 2 -d 0.2 -p $process_id`  #`ps -o psr $process_id`
const data_filename_ram  = "Resources_Monitor/ram_"*file_name
const data_filename_cpu = "Resources_Monitor/cpu_"*file_name
const data_filename_all = "Resources_Monitor/all_"*file_name
const max_seconds   = max_hours*60.0*60.0
const update_time   = 1

time_stamp    = 0.0
ram_array     = Array{Float64}(undef,0)
cpu_array     = Array{Float64}(undef,0)
time_array    = Array{Float64}(undef,0)
ram_final     = 0.0
cpu_final     = 0.0


function check_string_number(a)
    tryparse(Float64, a) !== nothing
end

function h5write_save(file_name,ram_array,cpu_array,time_array)
    file_h5=h5open(file_name*".h5", "w")
    file_h5["ram"]  = ram_array
    file_h5["cpu"]  = cpu_array
    file_h5["time"] = time_array
    close(file_h5)
end

println("#Starting monitoring RAM")
while time_stamp<max_seconds
    inp=Pipe()
    out=Pipe()
    err=Pipe()
    inp2=Pipe()
    out2=Pipe()
    err2=Pipe()
    try
        run(pipeline(pipeline(cmd_ram, `tail -n 1`), stdin=inp, stdout=out, stderr=err), wait=false)
        run(pipeline(pipeline(cmd_cpu, `tail -n 1`, `awk '{print $9}'`), stdin=inp2, stdout=out2, stderr=err2), wait=false)
    catch error_catch
        close(out.in)
        close(err.in)
        close(inp)
        close(out2.in)
        close(err2.in)
        close(inp2)
        println("Exiting the code: ", error_catch)
        error("Finished")#exit()
    end
    close(out.in)
    close(err.in)
    string_out = String(read(out))
    close(inp)
    #
    close(out2.in)
    close(err2.in)
    string_out2 = String(read(out2))
    close(inp2)
    #
    global time_stamp = time() - time_start
    global time_array = vcat(time_array, time_stamp)
    #
    string_temp_ram=string_out
    string_temp_cpu=string_out2
    length(string_temp_ram)==0 || length(string_temp_ram)==0 ? error("Finished. Time: ", time_stamp) : nothing
    #
    try
        global ram_final = parse(Float64, string_temp_ram)/1024^2
        global cpu_final = parse(Float64, string_temp_cpu)/100
        global ram_array  = vcat(ram_array,  ram_final )
        global cpu_array  = vcat(cpu_array,  cpu_final )
    catch
        println("#Finished monitoring")
        break
    end
    #
    println(ram_final,"  ",cpu_final,"  ",time_stamp)
    #
    #h5write_save(data_filename_all, ram_array,cpu_array,time_array)
    #
    my_plot  = plot(time_array, ram_array  , legend = false, xlabel="Time (s)", ylabel="RAM (GB)")
    my_plot2 = plot(time_array,cpu_array, legend = false, xlabel="Time (s)", ylabel="Cores")
    png(my_plot, data_filename_ram)
    png(my_plot2, data_filename_cpu)
    sleep(update_time)
end

error("#Run out of time")
