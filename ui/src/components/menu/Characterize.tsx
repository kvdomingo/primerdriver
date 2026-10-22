import { useForm, useSelector } from "@tanstack/react-form";
import { useMutation } from "@tanstack/react-query";
import { MDBCol as Col, MDBRow as Row } from "mdbreact";
import api from "../../api";
import {
  DnaSequenceInput,
  Form,
  MutationType,
  NumberMismatch,
  characterizeSchema,
} from "../form";
import Result from "./Result";

function Characterize() {
  const mode = "CHAR" as const;
  const primerDriver = useMutation({
    mutationFn: async (body: object) => (await api.data.primerDriver(body)).data,
  });

  const form = useForm({
    defaultValues: { sequence: "", mismatched_bases: 0, mutation_type: "" },
    validators: { onSubmit: characterizeSchema },
    onSubmit: ({ value }) => primerDriver.mutate({ mode, ...value }),
  });
  const isValid = useSelector(
    form.store,
    state => characterizeSchema.safeParse(state.values).success,
  );

  return (
    <>
      <Form
        title="Primer Characterization"
        onSubmit={form.handleSubmit}
        isValid={isValid}
      >
        <form.Field name="sequence">
          {field => (
            <DnaSequenceInput
              sequence={field.state.value}
              onChange={value =>
                field.handleChange(value.toUpperCase().replace(/[^ACGT]/g, ""))
              }
            />
          )}
        </form.Field>
        <Row className="row-cols-1 row-cols-md-2">
          <Col>
            <form.Field name="mismatched_bases">
              {field => (
                <NumberMismatch
                  handleChangeInt={field.handleChange}
                  mismatched_bases={field.state.value}
                />
              )}
            </form.Field>
          </Col>
          <Col>
            <form.Field name="mutation_type">
              {field => (
                <MutationType
                  handleChange={field.handleChange}
                  mutation_type={field.state.value}
                />
              )}
            </form.Field>
          </Col>
        </Row>
      </Form>
      <Result mode={mode} status={primerDriver.status} data={primerDriver.data} />
    </>
  );
}

export default Characterize;
