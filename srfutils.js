#!/usr/bin/env k8

/*******************************
 * Command line option parsing *
 *******************************/

var getopt = function(args, ostr) {
	var oli; // option letter list index
	if (typeof(getopt.place) == 'undefined')
		getopt.ind = 0, getopt.arg = null, getopt.place = -1;
	if (getopt.place == -1) { // update scanning pointer
		if (getopt.ind >= args.length || args[getopt.ind].charAt(getopt.place = 0) != '-') {
			getopt.place = -1;
			return null;
		}
		if (getopt.place + 1 < args[getopt.ind].length && args[getopt.ind].charAt(++getopt.place) == '-') { // found "--"
			++getopt.ind;
			getopt.place = -1;
			return null;
		}
	}
	var optopt = args[getopt.ind].charAt(getopt.place++); // character checked for validity
	if (optopt == ':' || (oli = ostr.indexOf(optopt)) < 0) {
		if (optopt == '-') return null; //  if the user didn't specify '-' as an option, assume it means null.
		if (getopt.place < 0) ++getopt.ind;
		return '?';
	}
	if (oli+1 >= ostr.length || ostr.charAt(++oli) != ':') { // don't need argument
		getopt.arg = null;
		if (getopt.place < 0 || getopt.place >= args[getopt.ind].length) ++getopt.ind, getopt.place = -1;
	} else { // need an argument
		if (getopt.place >= 0 && getopt.place < args[getopt.ind].length)
			getopt.arg = args[getopt.ind].substr(getopt.place);
		else if (args.length <= ++getopt.ind) { // no arg
			getopt.place = -1;
			if (ostr.length > 0 && ostr.charAt(0) == ':') return ':';
			return '?';
		} else getopt.arg = args[getopt.ind]; // white space
		getopt.place = -1;
		++getopt.ind;
	}
	return optopt;
}

/************
 * Commands *
 ************/

function srf_cmd_enlong(args) {
	var c, opt = { min_len:200 };
	while ((c = getopt(args, "l:")) != null) {
		if (c == 'l') opt.min_len = parseInt(getopt.arg);
	}
	if (args.length == getopt.ind) {
		print("Usage: srfutils.js enlong [options] <circ.fa>");
		print("Options:");
		print("  -l INT    min length [" + opt.min_len + "]");
		return 1;
	}
	var buf = new Bytes();
	var file = new File(args[getopt.ind]);

	function enlong_seq(seq) {
		var t = seq;
		while (t.length < opt.min_len)
			t += seq;
		return t;
	}

	var name = null, comm = null, seq = null;
	while (file.readline(buf) >= 0) {
		var m, l = buf.toString();
		if ((m = /^>(\S+)(.*)/.exec(l)) != null) {
			if (name != null) {
				print(">" + name + comm);
				print(enlong_seq(seq));
			}
			name = m[1], comm = m[2], seq = "";
		} else {
			seq += l;
		}
	}
	if (name != null) {
		print(">" + name + comm);
		print(enlong_seq(seq));
	}
	file.close();
	buf.destroy();
}

function srf_cmd_filter(args) {
	var c, opt = { min_len:200, min_cov:0.9, min_iden:0.9 };
	while ((c = getopt(args, "l:c:d:")) != null) {
		if (c == 'l') opt.min_len = parseInt(getopt.arg);
		else if (c == 'c') opt.min_cov = parseFloat(getopt.arg);
		else if (c == 'd') opt.min_iden = parseFloat(getopt.arg);
	}
	if (args.length == getopt.ind) {
		print("Usage: srfutils.js filter [options] <in.paf>");
		print("Options:");
		print("  -l INT      min alignment length [" + opt.min_len + "]");
		print("  -c FLOAT    min coverage [" + opt.min_cov + "]");
		print("  -d FLOAT    min identity [" + opt.min_iden + "]");
		print("Notes: suggested minimap2 setting:");
		print("  minimap2 -c -N1000000 -f1000 -r100,100 <(srfutils.js enlong circ.fa) reads.fa");
		return 1;
	}
	var buf = new Bytes();
	var file = new File(args[getopt.ind]);
	while (file.readline(buf) >= 0) {
		var line = buf.toString();
		var t = line.split("\t", 12);
		for (var i = 1; i <= 3; ++i) t[i] = parseInt(t[i]);
		for (var i = 6; i < t.length; ++i) t[i] = parseInt(t[i]);
		if (t[8] - t[7] < opt.min_len) continue;
		var clip = [t[7], t[6] - t[8]];
		if (t[4] == "+") {
			clip[0] = clip[0] < t[2]? clip[0] : t[2];
			clip[1] = clip[1] < t[1] - t[3]? clip[1] : t[1] - t[3];
		} else {
			clip[1] = clip[1] < t[2]? clip[1] : t[2];
			clip[0] = clip[0] < t[1] - t[3]? clip[0] : t[1] - t[3];
		}
		if (t[8] - t[7] < (t[8] - t[7] + clip[0] + clip[1]) * opt.min_cov) continue;
		if (t[9] < t[10] * opt.min_iden) continue;
		print(line);
	}
	file.close();
	buf.destroy();
}

/*************************
 ***** main function *****
 *************************/

function main(args)
{
	if (args.length == 0) {
		print("Usage: srfutils.js <command> [arguments]");
		print("Commands:");
		print("  enlong    enlong circ contigs");
		print("  filter    filter read-to-circ alignment");
		exit(1);
	}

	var cmd = args.shift();
	if (cmd == 'filter') srf_cmd_filter(args);
	else if (cmd == 'enlong') srf_cmd_enlong(args);
	else throw Error("unrecognized command: " + cmd);
}

main(arguments);
