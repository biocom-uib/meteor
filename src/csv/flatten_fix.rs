use serde::{Serializer, ser::SerializeStruct};

pub trait SerializeFlat {
    const FIELD_COUNT: usize;

    fn serialize_flat<Ser: SerializeStruct>(&self, row: &mut Ser) -> Result<(), Ser::Error>;
}

pub fn serialize_flat_struct<Ser: Serializer, T: SerializeFlat>(
    serializer: Ser,
    constructor: &'static str,
    x: &T,
) -> Result<Ser::Ok, Ser::Error> {
    let mut row = serializer.serialize_struct(constructor, T::FIELD_COUNT)?;
    x.serialize_flat(&mut row)?;
    row.end()
}
