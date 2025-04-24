use std::{io::BufRead, marker::PhantomData};

use quick_xml::{Reader, events::{Event, BytesStart}};


type XmlResult<T> = quick_xml::Result<T>;

pub enum SubfieldsAction {
    Skip,
    Populate,
}

#[allow(unused_variables)]
pub trait EntryBuilder: Default {
    type Entry;

    // accession
    fn accession(&mut self, acc: &str) {}

    // name
    fn name(&mut self, name: &str) {}

    // sequence
    fn sequence(&mut self, sequence: &str) {}

    // dbreference
    #[must_use]
    fn begin_dbreference(&mut self, ref_type: &str, ref_id: &str) -> SubfieldsAction {
        SubfieldsAction::Skip
    }

    fn dbreference_property(&mut self, prop_type: &str, prop_value: &str) {}
    fn finish_dbreference(&mut self) {}

    // organism
    #[must_use]
    fn begin_organism(&mut self) -> SubfieldsAction {
        SubfieldsAction::Skip
    }

    fn organism_scientific_name(&mut self, name: &str) {}
    fn organism_common_name(&mut self, name: &str) {}

    //  organism/dbReference
    #[must_use]
    fn begin_organism_dbreference(&mut self, ref_type: &str, ref_id: &str) -> SubfieldsAction {
        SubfieldsAction::Skip
    }
    fn organism_dbreference_property(&mut self, prop_type: &str, prop_value: &str) {}
    fn finish_organism_dbreference(&mut self) {}

    fn finish_organism(&mut self) {}

    // finish building
    fn finish(self) -> Option<Self::Entry>;
}

#[derive(Default)]
pub struct IdentityBuilder<EB>(EB);

impl<EB: EntryBuilder> EntryBuilder for IdentityBuilder<EB> {
    type Entry = EB;

    fn accession(&mut self, acc: &str) {
        self.0.accession(acc)
    }

    fn name(&mut self, name: &str) {
        self.0.name(name)
    }

    fn sequence(&mut self, sequence: &str) {
        self.0.sequence(sequence)
    }

    fn begin_dbreference(&mut self, ref_type: &str, ref_id: &str) -> SubfieldsAction {
        self.0.begin_dbreference(ref_type, ref_id)
    }

    fn dbreference_property(&mut self, prop_type: &str, prop_value: &str) {
        self.0.dbreference_property(prop_type, prop_value)
    }

    fn finish_dbreference(&mut self) {
        self.0.finish_dbreference()
    }

    fn begin_organism(&mut self) -> SubfieldsAction {
        self.0.begin_organism()
    }

    fn organism_scientific_name(&mut self, name: &str) {
        self.0.organism_scientific_name(name)
    }

    fn organism_common_name(&mut self, name: &str) {
        self.0.organism_common_name(name)
    }

    fn begin_organism_dbreference(&mut self, ref_type: &str, ref_id: &str) -> SubfieldsAction {
        self.0.begin_organism_dbreference(ref_type, ref_id)
    }

    fn organism_dbreference_property(&mut self, prop_type: &str, prop_value: &str) {
        self.0.organism_dbreference_property(prop_type, prop_value)
    }

    fn finish_organism_dbreference(&mut self) {
        self.0.finish_organism_dbreference()
    }

    fn finish_organism(&mut self) {
        self.0.finish_organism()
    }

    // finish building
    fn finish(self) -> Option<Self::Entry> {
        Some(self.0)
    }
}


pub struct UniProtXmlReader<R: BufRead> {
    xml_reader: Reader<R>,
    buffer: Vec<u8>,
    skip_to_end_buffer: Vec<u8>,
}

pub mod dbref_types {
    pub const NCBI_TAXONOMY: &str = "NCBI Taxonomy";
}

#[allow(dead_code)]
mod attrs {
    pub const TYPE: &str  = "type";
    pub const ID: &str    = "id";
    pub const VALUE: &str = "value";

    pub mod u8 {
        pub const TYPE:  &[u8] = super::TYPE.as_bytes();
        pub const ID:    &[u8] = super::ID.as_bytes();
        pub const VALUE: &[u8] = super::VALUE.as_bytes();
    }
}

#[allow(dead_code)]
mod tags {
    pub const ACCESSION: &str   = "accession";
    pub const COPYRIGHT: &str   = "copyright";
    pub const COMMON: &str      = "common";
    pub const DBREFERENCE: &str = "dbReference";
    pub const ENTRY: &str       = "entry";
    pub const NAME: &str        = "name";
    pub const ORGANISM: &str    = "organism";
    pub const PROPERTY: &str    = "property";
    pub const SCIENTIFIC: &str  = "scientific";
    pub const SEQUENCE: &str    = "sequence";
    pub const UNIPROT: &str     = "uniprot";

    pub mod u8 {
        pub const ACCESSION: &[u8]   = super::ACCESSION.as_bytes();
        pub const COPYRIGHT: &[u8]   = super::COPYRIGHT.as_bytes();
        pub const COMMON: &[u8]      = super::COMMON.as_bytes();
        pub const DBREFERENCE: &[u8] = super::DBREFERENCE.as_bytes();
        pub const ENTRY: &[u8]       = super::ENTRY.as_bytes();
        pub const NAME: &[u8]        = super::NAME.as_bytes();
        pub const ORGANISM: &[u8]    = super::ORGANISM.as_bytes();
        pub const PROPERTY: &[u8]    = super::PROPERTY.as_bytes();
        pub const SCIENTIFIC: &[u8]  = super::SCIENTIFIC.as_bytes();
        pub const SEQUENCE: &[u8]    = super::SEQUENCE.as_bytes();
        pub const UNIPROT: &[u8]     = super::UNIPROT.as_bytes();
    }
}

macro_rules! set_from_read_text {
    ($self:ident, $builder:ident, $tag:ident, $setter:ident) => {{
        $self.read_text(tags::u8::$tag, |text| {
            $builder.$setter(text);
            Ok(())
        })?;
    }}
}

macro_rules! expect_end_of {
    ($start:ident, $end:ident) => {{
        assert!(tags::u8::$start == $end.name().as_ref());
    }}
}

macro_rules! skip_to_end {
    ($self:ident, $start:ident) => {{
        $self.xml_reader.read_to_end_into($start.to_end().name(), &mut $self.skip_to_end_buffer)?;
        $self.skip_to_end_buffer.clear();
    }}
}

macro_rules! read_dbreference {
    ($self:ident, $builder:ident, $start:ident, $begin:ident, $add_prop:ident, $finish:ident) => {{
        let action = UniProtXmlReader::begin_read_dbreference(&$self.xml_reader, &$start, |ref_type, ref_id| {
            $builder.$begin(ref_type, ref_id)
        })?;

        if let Some(SubfieldsAction::Populate) = action {
            $self.read_dbreference_properties(|prop_type, prop_value| {
                $builder.$add_prop(prop_type, prop_value)
            })?;
        } else {
            skip_to_end!($self, $start);
        }

        $builder.$finish();
    }}
}


impl<R: BufRead> UniProtXmlReader<R> {
    pub fn new(reader: R) -> Self {
        let mut xml_reader = Reader::from_reader(reader);

        xml_reader.config_mut().expand_empty_elements = true;
        xml_reader.config_mut().trim_text(true);

        Self {
            xml_reader,
            buffer: Vec::new(),
            skip_to_end_buffer: Vec::new(),
        }
    }

    fn read_text<OnText>(&mut self, tag: &[u8], mut on_text: OnText) -> XmlResult<()>
    where
        OnText: FnMut(&str) -> XmlResult<()>,
    {
        loop {
            match self.xml_reader.read_event_into(&mut self.buffer)? {
                Event::Start(e) => {
                    skip_to_end!(self, e);
                }
                Event::Text(e) => {
                    let unescaped = e.unescape()?;
                    on_text(&unescaped)?;
                }
                Event::End(e) => {
                    assert!(e.name().as_ref() == tag);
                    break;
                }
                Event::Eof => break,
                _ => {},
            }
        }

        Ok(())
    }

    fn find_next_entry_start(&mut self) -> XmlResult<bool> {
        loop {
            match self.xml_reader.read_event_into(&mut self.buffer)? {
                Event::DocType(_) => {}
                Event::Decl(_) => {}
                Event::Start(e) => match e.name().as_ref() {
                    tags::u8::ENTRY => return Ok(true),
                    tags::u8::UNIPROT => {
                        // continue (enter)
                    }
                    tags::u8::COPYRIGHT => {
                        skip_to_end!(self, e);
                    }
                    name => {
                        eprintln!("Unknown tag found while scanning entries: {:?}", std::str::from_utf8(name));
                        skip_to_end!(self, e);
                    }
                },
                Event::End(e) => {
                    expect_end_of!(UNIPROT, e);
                }
                Event::Eof => {
                    return Ok(false);
                }
                evt => {
                    eprintln!("Unexpected event while scanning entries: {evt:?}");
                }
            };

            self.buffer.clear();
        }
    }

    pub fn into_entries<EB>(self) -> EntriesIter<R, EB> {
        EntriesIter {
            reader: self,
            error: false,
            marker: PhantomData,
        }
    }

    fn read_next_entry<EB: EntryBuilder>(&mut self) -> XmlResult<Option<EB::Entry>> {
        loop {
            if self.find_next_entry_start()? {
                let mut builder = EB::default();
                self.read_entry(&mut builder)?;

                if let Some(entry) = builder.finish() {
                    return Ok(Some(entry))
                }
            } else {
                return Ok(None)
            }
        }
    }

    fn read_entry<EB: EntryBuilder>(&mut self, builder: &mut EB) -> XmlResult<()> {
        loop {
            match self.xml_reader.read_event_into(&mut self.buffer)? {
                Event::Start(e) => {
                    match e.name().as_ref() {
                        tags::u8::ACCESSION => {
                            set_from_read_text!(self, builder, ACCESSION, accession);
                        }
                        tags::u8::NAME => {
                            set_from_read_text!(self, builder, NAME, name);
                        }
                        tags::u8::DBREFERENCE => {
                            read_dbreference!(self, builder, e,
                                begin_dbreference,
                                dbreference_property,
                                finish_dbreference
                            );
                        }
                        tags::u8::SEQUENCE => {
                            set_from_read_text!(self, builder, SEQUENCE, sequence);
                        }
                        tags::u8::ORGANISM => {
                            if let SubfieldsAction::Populate = builder.begin_organism() {
                                self.read_organism(builder)?;
                            } else {
                                skip_to_end!(self, e);
                            }
                            builder.finish_organism();
                        }
                        _ => {
                            skip_to_end!(self, e);
                        }
                    };
                }
                Event::End(e) => {
                    expect_end_of!(ENTRY, e);
                    break;
                }
                evt => {
                    eprintln!("Unexpected event while reading <entry>: {evt:?}");
                }
            }

            self.buffer.clear();
        }

        Ok(())
    }

    fn read_organism<EB: EntryBuilder>(&mut self, builder: &mut EB) -> XmlResult<()> {
        loop {
            match self.xml_reader.read_event_into(&mut self.buffer)? {
                Event::Start(e) => match e.name().as_ref() {
                    tags::u8::SCIENTIFIC => {
                        set_from_read_text!(self, builder, SCIENTIFIC, organism_scientific_name);
                    }
                    tags::u8::COMMON => {
                        set_from_read_text!(self, builder, NAME, organism_common_name);
                    }
                    tags::u8::DBREFERENCE => {
                        read_dbreference!(self, builder, e,
                            begin_organism_dbreference,
                            organism_dbreference_property,
                            finish_organism_dbreference
                        );
                    }
                    _ => {
                        skip_to_end!(self, e);
                    }
                },
                Event::End(e) => {
                    expect_end_of!(ORGANISM, e);
                    break;
                }
                evt => {
                    eprintln!("Unexpected event while reading <organism>: {evt:?}",);
                }
            }

            self.buffer.clear();
        }

        Ok(())
    }

    fn begin_read_dbreference(
        xml_reader: &Reader<R>,
        start: &BytesStart<'_>,
        begin: impl FnOnce(&str, &str) -> SubfieldsAction,
    ) -> XmlResult<Option<SubfieldsAction>> {

        let decoder = xml_reader.decoder();

        let mut ref_type = None;
        let mut id = None;

        for attr in start.attributes() {
            let attr = attr?;

            match attr.key.as_ref() {
                b"type" => ref_type = Some(attr.value),
                b"id" => id = Some(attr.value),
                _ => {}
            };
        }

        if let (Some(ref_type), Some(id)) = (ref_type, id) {
            let ref_type = decoder.decode(&ref_type)?;
            let id = decoder.decode(&id)?;
            Ok(Some(begin(&ref_type, &id)))
        } else {
            Ok(None)
        }
    }

    fn read_dbreference_properties(
        &mut self,
        mut add_prop: impl FnMut(&str, &str),
    ) -> XmlResult<()> {
        let decoder = self.xml_reader.decoder();

        loop {
            #[allow(clippy::single_match)]
            match self.xml_reader.read_event_into(&mut self.buffer)? {
                Event::Start(e) => {
                    match e.name().as_ref() {
                        tags::u8::PROPERTY => {
                            let mut prop_type = None;
                            let mut prop_value = None;

                            for attr in e.attributes() {
                                let attr = attr?;

                                match attr.key.as_ref() {
                                    attrs::u8::TYPE => prop_type = Some(attr.value),
                                    attrs::u8::VALUE => prop_value = Some(attr.value),
                                    _ => {}
                                };
                            }

                            if let (Some(prop_type), Some(prop_value)) = (prop_type, prop_value) {
                                let prop_type = decoder.decode(&prop_type)?;
                                let prop_value = decoder.decode(&prop_value)?;
                                add_prop(&prop_type, &prop_value);
                            }
                        }
                        _ => {}
                    }

                    skip_to_end!(self, e);
                }
                Event::End(e) => {
                    expect_end_of!(DBREFERENCE, e);
                    break;
                }
                evt => {
                    eprintln!("Unexpected event while reading <dbReference>: {evt:?}");
                }
            }

            self.buffer.clear();
        }

        Ok(())
    }
}

pub struct EntriesIter<R: BufRead, EB> {
    reader: UniProtXmlReader<R>,
    error: bool,
    marker: PhantomData<fn() -> EB>,
}

impl<R, EB> Iterator for EntriesIter<R, EB>
where
    R: BufRead,
    EB: EntryBuilder,
{
    type Item = XmlResult<EB::Entry>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.error {
            return None;
        }

        let r = self.reader.read_next_entry::<EB>();

        if r.is_err() {
            self.error = true;
        }

        r.transpose()
    }
}
