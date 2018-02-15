PreFLAGS = -classpath $(CURDIR)
PostFLAGS = -Xlint:unchecked
JC = javac
.SUFFIXES: .java .class
.java.class:
	$(JC) $(PreFLAGS) $*.java $(PostFLAGS)

MESHI = \
        $.meshi/*.java \
	$.meshi/*/*.java \
	$.meshi/*/*/*.java \
	$.meshi/*/*/*/*.java \
	$.meshi/*/*/*/*/*.java \
	$.programs/*.java


default: classes

classes: $(MESHI:.java=.class)

clean:
	$(RM) *.class
	$(RM) */*.class
	$(RM) */*/*.class
	$(RM) */*/*/*.class
	$(RM) */*/*/*/*.class
	$(RM) */*/*/*/*/*.class
	$(RM) */*/*/*/*/*/*.class
